"""Suppress line numbering in footers of a .docx (Word does not number footer lines, but some
renderers do). Inserts <w:suppressLineNumbers/> at the schema-correct position of each footer pPr."""
import re
import shutil
import sys
import zipfile
from pathlib import Path

src = Path(sys.argv[1])
tmp = src.with_suffix(".tmp.docx")
with zipfile.ZipFile(src) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if re.match(r"word/footer\d*\.xml$", item.filename):
            x = data.decode("utf8")
            x = x.replace("<w:pPr><w:jc", "<w:pPr><w:suppressLineNumbers/><w:jc")
            data = x.encode("utf8")
        zout.writestr(item, data)
shutil.move(tmp, src)
print("post-processed", src)
