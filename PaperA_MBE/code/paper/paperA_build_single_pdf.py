"""Build MBE single-PDF submission via LibreOffice + fitz merge."""
import subprocess, os, io, shutil, sys

BASE = r"d:\人类正选择基因项目"
SOFFICE = r"C:\Program Files\LibreOffice\program\soffice.exe"
TMP = os.path.join(BASE, "results", "paper", "_single_pdf")
os.makedirs(TMP, exist_ok=True)
log = []

def conv(docx, outdir):
    r = subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf", "--outdir", outdir, docx],
                       capture_output=True, text=True, timeout=300)
    pdf = os.path.join(outdir, os.path.splitext(os.path.basename(docx))[0] + ".pdf")
    if not os.path.exists(pdf):
        raise SystemExit("convert failed: %s\n%s" % (docx, (r.stdout + r.stderr)[-800:]))
    log.append("converted %s -> %s (%d B)" % (os.path.basename(docx), os.path.basename(pdf), os.path.getsize(pdf)))
    return pdf

# MBE formatting pass on the manuscript docx (line numbers, double spacing, 25 mm margins)
subprocess.run([sys.executable, os.path.join(BASE, "scripts", "paperA_format_docx.py"),
                os.path.join(BASE, "results", "paper", "manuscript_english_v9.docx")], check=True)
log.append("docx MBE formatting applied (line numbers / double spacing / 25mm margins)")

ms_pdf = conv(os.path.join(BASE, "results", "paper", "manuscript_english_v9.docx"), TMP)
tb_pdf = conv(os.path.join(BASE, "results", "paper", "Main_Tables_v9.docx"), TMP)
cl_pdf = conv(os.path.join(BASE, "results", "paper", "cover_letter_v1.docx"), TMP)
shutil.copy2(cl_pdf, os.path.join(BASE, "results", "paper", "cover_letter_v1.pdf"))
log.append("cover letter PDF copied to results/paper/")

import fitz
out_pdf = os.path.join(BASE, "MBE_submission_PaperA_single.pdf")
doc = fitz.open(ms_pdf)
log.append("manuscript pages: %d (line numbers rendered natively from docx w:lnNumType)" % doc.page_count)
doc.insert_pdf(fitz.open(tb_pdf))
for f in ["Figure1_analysis_pipeline.pdf", "Figure2_BUSTED_pvalue_distribution.pdf",
          "Figure3_GDS_RDS_scatter.pdf", "Figure4_LOO_enrichment_matrix.pdf",
          "Figure5_independent_evidence.pdf", "Figure6_class_functional_content.pdf",
          "Figure7_integrative_answer.pdf"]:
    fd = fitz.open(os.path.join(BASE, "MBE_submission_PaperA", "03_figures", f))
    doc.insert_pdf(fd)
    fd.close()
    log.append("appended " + f)
doc.set_metadata({"title": "Coding versus Regulatory Selection in Human Evolution",
                  "author": "Xianming Wu, Li Zhang"})
doc.save(out_pdf, deflate=True)
doc.close()
d = fitz.open(out_pdf)
log.append("FINAL: %s (%d pages, %.2f MB)" % (os.path.basename(out_pdf), d.page_count, os.path.getsize(out_pdf)/1e6))
d.close()
io.open(r"C:\Users\admin\.workbuddy\tmp_single.txt", "w", encoding="utf-8").write("\n".join(log))
print("done")
