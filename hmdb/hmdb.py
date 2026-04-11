import csv
from lxml import etree

def get_text_or_na(element, path, ns=None):
    """
    安全获取元素的文本内容（支持命名空间）
    :param element: 父元素
    :param path: XPath表达式
    :param ns: 命名空间字典
    :return: 文本内容或'NA'
    """
    node = element.find(path, namespaces=ns) if ns else element.find(path)
    return node.text if node is not None and node.text else 'NA'

def parse_large_hmdb(xml_file, output_csv=None):
    """
    高效解析大型HMDB XML文件
    :param xml_file: 输入XML文件路径
    :param output_csv: 输出CSV文件路径（None则返回生成器）
    :return: 生成器或None（直接写入文件）
    """
    ns = {'hmdb': 'http://www.hmdb.ca'}
    fieldnames = [
        'accession', 'monisotopic_molecular_weight', 'iupac_name', 'name',
        'chemical_formula', 'cas_registry_number', 'smiles', 'kingdom',
        'direct_parent', 'super_class', 'class', 'sub_class',
        'molecular_framework'
    ]

    # 直接写入CSV模式
    if output_csv:
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for event, elem in etree.iterparse(xml_file, tag='{http://www.hmdb.ca}metabolite'):
                row = {
                    'accession': get_text_or_na(elem, 'hmdb:accession', ns),
                    'monisotopic_molecular_weight': get_text_or_na(elem, 'hmdb:monisotopic_molecular_weight', ns),
                    'iupac_name': get_text_or_na(elem, 'hmdb:iupac_name', ns),
                    'name': get_text_or_na(elem, 'hmdb:name', ns),
                    'chemical_formula': get_text_or_na(elem, 'hmdb:chemical_formula', ns),
                    'cas_registry_number': get_text_or_na(elem, 'hmdb:cas_registry_number', ns),
                    'smiles': get_text_or_na(elem, 'hmdb:smiles', ns),
                    'kingdom': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:kingdom', ns),
                    'direct_parent': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:direct_parent', ns),
                    'super_class': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:super_class', ns),
                    'class': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:class', ns),
                    'sub_class': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:sub_class', ns),
                    'molecular_framework': get_text_or_na(elem, 'hmdb:taxonomy/hmdb:molecular_framework', ns),
                }
                writer.writerow(row)
                
                # 内存清理
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        print(f"数据已写入 {output_csv}")
        
        
print("开始解析HMDB XML文件...")
parse_large_hmdb('hmdb_metabolites.xml', output_csv='hmdb_output.csv')
print("解析完成！结果已保存到 hmdb_output.csv")