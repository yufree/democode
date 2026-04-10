"""
litdigest.py — 文献追踪 + LLM 评分周报生成器

用法:
    python litdigest.py

配置:
    - RSS_FEEDS: 订阅源列表（支持 PubMed / ACS 等）
    - DAYS: 抓取最近几天的文章
    - TOP_N: 最终输出 Top N 篇
    - LLM_BACKEND: "ollama" 或 "openai"

API Key（OpenAI 模式）:
    export OPENAI_API_KEY=sk-...
"""

import os
import re
import time
import datetime
import logging
import feedparser
import requests
from bs4 import BeautifulSoup

# ─── 配置区 ────────────────────────────────────────────────────────────────────

RSS_FEEDS = [
    # PubMed 搜索订阅（摘要直接在 RSS 里，无需爬页面）
    "https://pubmed.ncbi.nlm.nih.gov/rss/search/1roWqUHnCMOw0LYZD2R1suracmWnrHoePMKeDCnAP7yHFdPILE/?limit=20&utm_campaign=pubmed-2&fc=20240425120030",
    # ACS 期刊（ES&T）
    "https://pubs.acs.org/action/showFeed?type=axatoc&feed=rss&jc=esthag",
]

DAYS = 7          # 抓取最近几天
TOP_N = 10        # 最终输出 Top N
SCORE_WEIGHTS = {"research": 0.618, "impact": 0.382}  # 黄金比例权重

# LLM 后端: "ollama" 使用本地模型，"openai" 使用 OpenAI API
LLM_BACKEND = "ollama"
OLLAMA_URL = "http://localhost:11434"
OLLAMA_MODEL = "deepseek-r1:8b"      # 本地模型名
OPENAI_MODEL = "gpt-4o-mini"

# ─── 日志 ──────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── RSS 抓取 ──────────────────────────────────────────────────────────────────

def fetch_recent_papers(feeds: list[str], days: int) -> list[dict]:
    """从 RSS 订阅源抓取最近 N 天的文章，返回标准化列表。"""
    cutoff = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(days=days)
    papers = []

    for url in feeds:
        log.info(f"Fetching: {url[:60]}...")
        try:
            feed = feedparser.parse(requests.get(url, timeout=10).text)
        except Exception as e:
            log.warning(f"Failed to fetch {url}: {e}")
            continue

        for entry in feed.entries:
            pub_date = _parse_date(entry)
            if pub_date and pub_date < cutoff:
                continue  # 跳过旧文章

            abstract = _extract_abstract(entry)
            papers.append({
                "title": entry.get("title", "No Title").strip(),
                "doi": entry.get("link") or entry.get("dc_identifier", ""),
                "date": pub_date.strftime("%Y-%m-%d") if pub_date else "unknown",
                "abstract": abstract,
            })

    log.info(f"Found {len(papers)} papers in the last {days} days.")
    return papers


def _parse_date(entry) -> datetime.datetime | None:
    for fmt in ("%a, %d %b %Y %H:%M:%S %z", "%a, %d %b %Y %H:%M:%S %Z"):
        raw = entry.get("published", "")
        try:
            dt = datetime.datetime.strptime(raw, fmt)
            return dt.astimezone(datetime.timezone.utc)
        except ValueError:
            pass
    if hasattr(entry, "published_parsed") and entry.published_parsed:
        try:
            return datetime.datetime(*entry.published_parsed[:6],
                                     tzinfo=datetime.timezone.utc)
        except Exception:
            pass
    return None


def _extract_abstract(entry) -> str:
    """优先从 RSS 条目取摘要，其次爬 DOI 页面。"""
    # PubMed: 摘要在 content 字段
    if hasattr(entry, "content") and entry.content:
        text = BeautifulSoup(entry.content[0].value, "html.parser").get_text()
        if len(text) > 100:
            return text.strip()

    # summary 字段
    summary = entry.get("summary", "")
    if len(summary) > 100:
        return BeautifulSoup(summary, "html.parser").get_text().strip()

    # 最后尝试爬 DOI 页面（ACS 等）
    doi = entry.get("link", "")
    if doi:
        return _scrape_abstract(doi)

    return ""


def _scrape_abstract(url: str) -> str:
    """爬取 ACS 期刊页面摘要，失败返回空字符串。"""
    try:
        time.sleep(2)
        resp = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=15)
        soup = BeautifulSoup(resp.content, "html.parser")
        # ACS 摘要选择器
        div = soup.find("div", class_="article_abstract-content", id="abstractBox")
        if div:
            p = div.find("p", class_="articleBody_abstractText")
            if p:
                return p.get_text(strip=True)
    except Exception as e:
        log.debug(f"Scrape failed for {url}: {e}")
    return ""

# ─── LLM 调用 ─────────────────────────────────────────────────────────────────

def score_paper(paper: dict) -> dict:
    """用 LLM 对单篇文章打分，返回含 research_score / impact_score 的字典。"""
    text = f"Title: {paper['title']}\n\nAbstract: {paper['abstract'][:1000]}"
    prompt = (
        "You are an environmental science expert. Score this paper on two dimensions:\n"
        "1. Research Score (0-100): innovation, rigor, reliability\n"
        "2. Social Impact Score (0-100): policy relevance, public attention, societal impact\n\n"
        f"{text}\n\n"
        "Reply ONLY in this format:\nResearch Score: <number>\nSocial Impact Score: <number>"
    )

    response = _call_llm(prompt)
    research = _extract_score(response, "Research Score")
    impact = _extract_score(response, "Social Impact Score")

    return {**paper, "research_score": research, "impact_score": impact,
            "total_score": research * SCORE_WEIGHTS["research"] + impact * SCORE_WEIGHTS["impact"]}


def _call_llm(prompt: str) -> str:
    if LLM_BACKEND == "openai":
        return _call_openai(prompt)
    return _call_ollama(prompt)


def _call_ollama(prompt: str) -> str:
    try:
        resp = requests.post(
            f"{OLLAMA_URL}/api/generate",
            json={"model": OLLAMA_MODEL, "prompt": prompt, "stream": False},
            timeout=60,
        )
        return resp.json().get("response", "")
    except Exception as e:
        log.error(f"Ollama error: {e}")
        return ""


def _call_openai(prompt: str) -> str:
    try:
        import openai
        client = openai.OpenAI(api_key=os.environ["OPENAI_API_KEY"])
        resp = client.chat.completions.create(
            model=OPENAI_MODEL,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=100,
        )
        return resp.choices[0].message.content.strip()
    except Exception as e:
        log.error(f"OpenAI error: {e}")
        return ""


def _extract_score(text: str, label: str) -> float:
    match = re.search(rf"{label}:\s*(\d+(?:\.\d+)?)", text)
    return float(match.group(1)) if match else 50.0  # 默认 50 分

# ─── 报告生成 ──────────────────────────────────────────────────────────────────

def generate_report(papers: list[dict], top_n: int) -> str:
    """生成 Markdown 格式的周报。"""
    ranked = sorted(papers, key=lambda x: x.get("total_score", 0), reverse=True)[:top_n]
    today = datetime.date.today().strftime("%Y-%m-%d")

    lines = [
        f"# 文献周报 {today}",
        f"\n> 来源：{len(RSS_FEEDS)} 个订阅源，最近 {DAYS} 天，共 {len(papers)} 篇，精选 Top {top_n}\n",
    ]

    for i, p in enumerate(ranked, 1):
        lines += [
            f"## {i}. {p['title']}",
            f"- **日期**: {p['date']}",
            f"- **DOI**: {p['doi']}",
            f"- **研究分**: {p['research_score']:.0f}  |  **影响分**: {p['impact_score']:.0f}  |  **综合分**: {p['total_score']:.1f}",
            f"\n{p['abstract'][:300]}...\n" if p['abstract'] else "",
        ]

    return "\n".join(lines)

# ─── 主流程 ───────────────────────────────────────────────────────────────────

def main():
    papers = fetch_recent_papers(RSS_FEEDS, DAYS)

    if not papers:
        log.warning("No papers found.")
        return

    log.info(f"Scoring {len(papers)} papers with {LLM_BACKEND}...")
    scored = []
    for i, p in enumerate(papers, 1):
        log.info(f"  [{i}/{len(papers)}] {p['title'][:50]}...")
        if p["abstract"]:
            scored.append(score_paper(p))
        else:
            log.warning(f"  Skipped (no abstract): {p['title'][:50]}")

    report = generate_report(scored, TOP_N)

    outfile = f"digest_{datetime.date.today()}.md"
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(report)

    log.info(f"Report saved: {outfile}")
    print(report)


if __name__ == "__main__":
    main()
