# from io import StringIO
from lxml import etree
import csv
xml = 'data/hmdb.xml'

context = etree.iterparse(xml, tag='metabolite')

csvfile = open('hmdb.csv', 'w')
fieldnames = ['accession', 'monisotopic_molecular_weight', 'iupac_name', 'name', 'chemical_formula', 'cas_registry_number', 'smiles', 'kingdom', 'direct_parent', 'super_class', 'class', 'sub_class', 'molecular_framework']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for event, elem in context:

    accession = elem.xpath('accession/text()')[0]
    try:
        monisotopic_molecular_weight = elem.xpath('monisotopic_molecular_weight/text()')[0]
    except:
        monisotopic_molecular_weight = 'NA'
    try:
        iupac_name = elem.xpath('iupac_name/text()')[0].encode('utf-8')
    except:
        iupac_name = 'NA'
    name = elem.xpath('name/text()')[0].encode('utf-8')
    try:
        chemical_formula = elem.xpath('chemical_formula/text()')[0]
    except:
        chemical_formula = 'NA'

    try:
        cas_registry_number = elem.xpath('cas_registry_number/text()')[0]
    except:
        cas_registry_number = 'NA'
    try:
        smiles = elem.xpath('smiles/text()')[0]
    except:
        smiles = 'NA'
    try:
        kingdom = elem.xpath('taxonomy/kingdom/text()')[0]
    except:
        kingdom = 'NA'
    try:
        direct_parent = elem.xpath('taxonomy/direct_parent/text()')[0]
    except:
        direct_parent = 'NA'
    try:
        super_class = elem.xpath('taxonomy/super_class/text()')[0]
    except:
        super_class = 'NA'
    try:
        classorg = elem.xpath('taxonomy/class/text()')[0]
    except:
        classorg = 'NA'
    try:
        sub_class = elem.xpath('taxonomy/sub_class/text()')[0]
    except:
        sub_class = 'NA'
    try:
        molecular_framework = elem.xpath('taxonomy/molecular_framework/text()')[0]
    except:
        molecular_framework = 'NA'

    writer.writerow({'accession': accession, 'monisotopic_molecular_weight': monisotopic_molecular_weight, 'iupac_name': iupac_name, 'name': name, 'chemical_formula': chemical_formula, 'cas_registry_number': cas_registry_number, 'smiles': smiles, 'kingdom': kingdom, 'direct_parent': direct_parent, 'super_class': super_class, 'class': classorg, 'sub_class': sub_class, 'molecular_framework': molecular_framework})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
del context
