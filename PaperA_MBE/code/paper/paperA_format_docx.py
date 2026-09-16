# -*- coding: utf-8 -*-
"""
Post-process a pandoc-generated docx for MBE submission formatting:
  - continuous line numbers (w:lnNumType, recommended by MBE)
  - double line spacing via docDefaults pPrDefault (MBE encourages double spaced)
  - 25 mm margins on all sides (1417 twips)
Usage: python paperA_format_docx.py <path-to.docx>
Modifies the file in place (a .bak backup is written alongside).
"""
import io, os, re, shutil, sys, zipfile

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

def patch(docx_path):
    shutil.copy2(docx_path, docx_path + '.bak')
    zin = zipfile.ZipFile(docx_path, 'r')
    items = {n: zin.read(n) for n in zin.namelist()}
    zin.close()

    doc = items['word/document.xml'].decode('utf-8')

    # --- page size + 25 mm margins + line numbers inside sectPr ---
    sect = ('<w:pgSz w:w="12240" w:h="15840"/>'
            '<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" '
            'w:header="720" w:footer="720" w:gutter="0"/>'
            '<w:lnNumType w:countBy="1" w:start="0" w:distance="240" w:restart="continuous"/>')
    if '<w:pgMar' in doc:
        def fix_pgmar(m):
            tag = m.group(0)
            for attr in ['top', 'right', 'bottom', 'left']:
                if 'w:%s="' % attr in tag:
                    tag = re.sub(r'w:%s="-?\d+"' % attr, 'w:%s="1417"' % attr, tag)
                else:
                    tag = tag.replace('/>', ' w:%s="1417"/>' % attr)
            return tag
        doc = re.sub(r'<w:pgMar[^>]*/>', fix_pgmar, doc)
        if 'w:lnNumType' not in doc:
            doc = re.sub(r'(<w:pgMar[^>]*/>)',
                         r'\1<w:lnNumType w:countBy="1" w:start="0" w:distance="240" w:restart="continuous"/>',
                         doc)
    elif '</w:footnotePr>' in doc:
        doc = doc.replace('</w:footnotePr>', '</w:footnotePr>\n    ' + sect, 1)
    elif '<w:sectPr>' in doc:
        doc = doc.replace('<w:sectPr>', '<w:sectPr>\n    ' + sect, 1)

    items['word/document.xml'] = doc.encode('utf-8')

    # --- double spacing via styles.xml docDefaults ---
    st = items['word/styles.xml'].decode('utf-8')
    if 'w:line="480"' not in st:
        if '<w:pPrDefault>' in st:
            st = re.sub(r'(<w:pPrDefault>\s*<w:pPr>)',
                        r'\1<w:spacing w:line="480" w:lineRule="auto"/>', st, count=1)
        elif '<w:docDefaults>' in st:
            st = st.replace('<w:docDefaults>',
                            '<w:docDefaults><w:pPrDefault><w:pPr><w:spacing w:line="480" w:lineRule="auto"/></w:pPr></w:pPrDefault>',
                            1)
    items['word/styles.xml'] = st.encode('utf-8')

    tmp = docx_path + '.tmp'
    zout = zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED)
    for n, data in items.items():
        zout.writestr(n, data)
    zout.close()
    os.replace(tmp, docx_path)
    print('formatted:', docx_path)

if __name__ == '__main__':
    patch(sys.argv[1])
