"""Build inst/registry/methods.csv and presets.csv.

Measured values (runtime, number of outputs, xCell2 preset results) come from
data-raw/smoke_results2.csv, written by validation/03_smoke_test.R. The curated
columns (licence, status, install specs, notes) are verified facts recorded
below, each with the evidence that supports it (see validation/smoke-test-*.md).

    python data-raw/build_registry.py
"""
import csv
from pathlib import Path

here = Path(__file__).resolve().parent
out_dir = here.parent / "inst" / "registry"

rows = list(csv.DictReader(open(here / "smoke_results2.csv", encoding="utf-8")))
multi = {r["method"]: r for r in rows if r["scenario"] == "multi"}
dup = {r["method"]: r for r in rows if r["scenario"] == "dup"}

VALIDATED = "R 4.5.3; Bioconductor 3.22; Windows 11; 2026-11-11"

CLASS_NOTE = {
    "fraction": "estimated composition; comparable across samples for the same cell type and reference",
    "abundance_score": "arbitrary units; comparable across samples for one cell type, not across cell types",
    "enrichment_score": "relative score within the run; not comparable across cell types or across studies",
    "abundance_estimate": "model coefficients proportional to cell abundance; comparable across samples for one cell type",
    "relative_quantity": "relative cell quantity; designed for expression changes",
    "fraction_modules": "data-dependent marker modules; map to cell types before comparing",
}

GH = "github::"
DATA_SHA = "e625e6c28ed14a30f9f40f159925cc9f0df4fa49"  # immunedeconv commit the reference data come from (data only)
# install spec syntax: type::name[@ref][#subdir];  type is cran, bioc or github
M = {
    "quantiseq": dict(n_out="11", secs="0.4", engine="native", species="human", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="fraction", packages="quadprog", install="",
        source="native port of quanTIseq (Finotello et al. 2019) after quantiseqr 1.18.0, with the TIL10 files and the HGNC table", version="bulkdeconv", pin="data:quantiseqr 1.18.0", license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (quantiseqr); notice in inst/COPYRIGHTS",
        expr_ev="tpm (run_quantiseq, is_arraydata=FALSE)", registration="not required",
        notes="10 immune cell types plus 'Other'; is_tumordata=TRUE | native: no ggplot2, tidyr, SummarizedExperiment, preprocessCore, limSolve or MASS needed (only the small package quadprog for the constrained least squares); the original looks gene symbols up in a column called Synonyms that its HGNC table does not have, so only previous symbols are used (same result)"),
    "mcp_counter": dict(n_out="10", secs="0.3", engine="native", install="", species="human", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="abundance_score", packages="",
        source="native port of MCP-counter (Becht et al. 2016); marker table from ebecht/MCPcounter@b6eac73", version="bulkdeconv", pin="data:MCPcounter@b6eac73",
        license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (MCPcounter); the marker table is bundled, so nothing is downloaded at run time", expr_ev="not stated; ran with TPM", registration="not required",
        notes="MCPcounter.estimate() downloads genes.txt at call time; the adapter fetches it once at the pinned commit and caches it | native: no curl or MCPcounter package, no download at run time (the original fetches genes.txt from GitHub on every call)"),
    "xcell": dict(n_out="67", secs="1.7", engine="native", install="", species="human", status="conditional", min_n="2", single="unsupported", indep="no (diff 0.566 when run with a duplicate)", cohort="yes",
        needs="", options="", cls="enrichment_score", packages="",
        source="native port of xCell (Aran et al. 2017) with its bundled signatures and spill-over matrix", version="bulkdeconv", pin="data:xCell@20e2919",
        license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (xCell); notice in inst/COPYRIGHTS", expr_ev="rnaseq=TRUE; ran with TPM", registration="not required",
        notes="scores depend on the other samples in the run; superseded by xcell2 | native: no GSVA, GSEABase, pracma, quadprog or digest needed; ssGSEA and NNLS are written out and checked"),
    "xcell2": dict(engine="native", species="human|mouse", status="ready", min_n="1", single="duplicate", indep="no (subtracts the cohort minimum)", cohort="yes",
        needs="reference preset", options="preset", cls="enrichment_score", packages="", install="",
        source="native xCell2 1.2.3 inference with its pre-trained reference presets", version="bulkdeconv", pin="data:xCell2 1.2.3 (38ae014)", license="GPL-3 (bulkdeconv)", lic_note="xCell2 GPL (>= 3); see inst/COPYRIGHTS for reference provenance",
        expr_ev="references trained on TPM or RMA; ran with TPM", registration="not required",
        notes="9 presets; at least 50% reference-gene overlap required. One sample is duplicated for compatibility; both copies then have the cohort minimum, so the scores are zero. Preset timings refer to the original implementation."),
    "estimate": dict(n_out="3", secs="0.4", engine="native", species="human", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="enrichment_score", packages="", install="",
        source="native port of ESTIMATE (Yoshihara et al. 2013) after tidyestimate 1.1.1; data from tidyestimate", version="bulkdeconv", pin="data:tidyestimate 1.1.1", license="GPL (>= 2) (ported into bulkdeconv, GPL-3)", lic_note="tidyestimate is GPL (>= 2); the ESTIMATE common genes and signatures are bundled with it (inst/COPYRIGHTS)",
        expr_ev="normalized expression; is_affymetrix=FALSE for RNA-seq", registration="not required",
        notes="stromal, immune and combined scores only; no tumour purity for RNA-seq | native: 0.4 s, same output as tidyestimate (tests/testthat/test-native.R)"),
    "consensus_tme": dict(n_out="19", secs="0.5", engine="native", install="", species="human", status="conditional", min_n="1", single="native", indep="no (diff 0.254 alone vs in cohort)", cohort="yes",
        needs="cancer_type (32 TCGA types or 'Unfiltered')", options="cancer_type", cls="enrichment_score", packages="",
        source="native port of ConsensusTME (Jimenez-Sanchez et al. 2019) with its precomputed consensus gene sets", version="bulkdeconv", pin="data:ConsensusTME@6e14ba3",
        license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (ConsensusTME); notice in inst/COPYRIGHTS", expr_ev="not stated; ran with TPM", registration="not required",
        notes="ssGSEA is relative within the run, so a single-sample value is not comparable with a cohort value; cancer_type 'Unfiltered' needs no cancer type | native: no GSVA, GSEABase or singscore needed; ssgsea only (the other statistics of the original are not offered)"),
    "immucellai2": dict(n_out="53", secs="1.7", engine="native", install="", species="human", status="ready", min_n="1", single="native", indep="nearly (diff 0.068)", cohort="no",
        needs="tissue ('tumor' by default)", options="tissue", cls="fraction", packages="",
        source="native port of ImmuCellAI 2.0 (Yang) from GuoBioinfoLab/ImmuCellAI-2.0@0135016; reference restricted to the 5510 marker genes", version="bulkdeconv", pin="data:ImmuCellAI-2.0@0135016",
        license="MIT (ImmuCellAI 2.0, ported)", lic_note="MIT + file LICENSE (Jingmin Yang, 2026); notice kept in inst/COPYRIGHTS", expr_ev="TPM, non-log, genes in rows (README)", registration="not required",
        notes="53 immune subtypes for tumour; 15-19 s on 4 samples | native: 1.7 s instead of 18 s (the original reads a 41 MB text reference on every call; here the reference is stored restricted to the 5510 marker genes, 2 MB)"),
    "imsig": dict(n_out="10", secs="0.3", engine="native", species="human", status="conditional", min_n="2", single="unsupported", indep="cohort-based by design (correlation filtering)", cohort="yes",
        needs="", options="", cls="abundance_score", packages="", install="",
        source="native port of ImSig (Nirmal et al. 2018) with its signature table", version="bulkdeconv", pin="data:imsig 1.1.3", license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (imsig); notice in inst/COPYRIGHTS",
        expr_ev="natural scale (documentation)", registration="not required",
        notes="needs 2 distinct samples (identical duplicates fail); outputs include non-cell signatures (Interferon, Proliferation, Translation); last release 2021-01 | native: no igraph, HiClimR, ggplot2, gridExtra or survival needed; the original sorts the samples by the T cell score, here they stay in input order"),
    "spatialdecon_safetme": dict(n_out="18", secs="1.1", engine="native", species="human", status="conditional", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="background matrix (constant 0.1 used)", options="", cls="abundance_estimate", packages="", install="",
        source="native port of SpatialDecon 1.20.1 (NanoString) with the bundled SAFE-TME matrix; logNormReg::lognlm written out", version="bulkdeconv", pin="data:SpatialDecon 1.20.1", license="MIT (SpatialDecon, ported)", lic_note="MIT (NanoString Technologies, 2021); notice kept in inst/COPYRIGHTS",
        expr_ev="normalized data plus a background matrix (GeoMx design)", registration="not required",
        notes="built for spatial data; bulk use is allowed by the package title but the meaning of the background for bulk is unvalidated; 'Unsuccessful convergence' warning on the test data; pulls in GeomxTools | native: no GeomxTools, SeuratObject or logNormReg needed; same coefficients as SpatialDecon; standard errors and p values are not computed"),
    "epic": dict(species="human", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="fraction", packages="EPIC",
        install=GH + "GfellerLab/EPIC@50a4f404f96c2842b2891b517b4e3bfaa6c64b8f",
        source="GitHub:GfellerLab/EPIC", version="1.1.7", pin="50a4f404f96c2842b2891b517b4e3bfaa6c64b8f",
        license="file LICENSE (academic, non-commercial)", lic_note="README: free for academic non-commercial use; for-profit use needs a licence",
        expr_ev="TPM or RPKM (README)", registration="not required",
        notes="no registration; a non-commercial notice is printed at run time; includes 'otherCells'"),
    "immucellai": dict(species="human", status="unavailable", min_n="", single="", indep="", cohort="",
        needs="", options="", cls="fraction", packages="", install="",
        source="GitHub:lydiaMyr/ImmuCellAI", version="0.1.0", pin="ad2877680cdd", license="GPL (version unstated)", lic_note="DESCRIPTION",
        expr_ev="FPKM or TPM (README)", registration="not required",
        notes="calls the removed GSVA API gsva(matrix, list, method=, parallel.sz=) and fails with GSVA 2.4.9 (Bioconductor 3.22); use immucellai2"),
    "mmcp_counter": dict(engine="native", install="", species="mouse", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="abundance_score", packages="",
        source="native port of mMCP-counter (Petitprez et al. 2020), genome GCRm39 signatures", version="bulkdeconv", pin="data:mMCP-counter@cec6971",
        license="GPL-3 (bulkdeconv)", lic_note="GPL-3 (mMCP-counter); notice in inst/COPYRIGHTS", expr_ev="not stated; ran with TPM", registration="not required",
        notes="signatures bundled for GCRm38 and GCRm39 (GCRm39 used); warning 'No genes were found for population(s): Basophils' on the test data | native: also works on one sample (the original fails because a one-column matrix loses its dimensions); only the GCRm39 symbols signature is bundled"),
    "dcq": dict(species="mouse", status="conditional", min_n="1", single="duplicate", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="relative_quantity", packages="ComICS", install="cran::ComICS",
        source="CRAN:ComICS", version="1.0.4", pin="cran-2026-11-11", license="GPL-2", lic_note="DESCRIPTION",
        expr_ev="designed for expression changes (delta); ran with TPM", registration="not required",
        notes="the bundled reference is a 61-gene ImmGen subset (immgen_dat, 207 populations) with 61 markers (59 present in the test data); 207 fine-grained outputs; last release 2018-05"),
    "immucellai_mouse": dict(n_out="37", secs="1.1", engine="native", install="", species="mouse", status="conditional", min_n="2", single="unsupported", indep="no (diff 2.62 when run with a duplicate)", cohort="yes",
        needs="", options="", cls="fraction", packages="",
        source="native port of ImmuCellAI-mouse 0.1.5 from BioInfoCloud/ImmuCellAI_mouse@b9fc2e7", version="bulkdeconv", pin="data:ImmuCellAI_mouse@b9fc2e7",
        license="MIT (ImmuCellAI-mouse, ported)", lic_note="MIT + file LICENSE (Yun Sen, 2025); notice kept in inst/COPYRIGHTS",
        expr_ev="FPKM or TPM for RNA-seq (README)", registration="not required",
        notes="37 outputs over three hierarchy layers; values depend strongly on the other samples | native: no GSVA or dplyr needed; ssGSEA and NNLS are written out and checked; one deliberate difference: the infiltration score matches sample names exactly instead of by regular expression (the original mixes samples when one name is part of another, e.g. S1 and S10)"),
    "ssmd": dict(species="mouse", status="conditional", min_n="6", single="unsupported", indep="cohort-based by design", cohort="yes",
        needs="tissue ('Inflammatory' uses the mouse cancer markers)", options="tissue", cls="fraction_modules", packages="SSMD|nnls|bcv",
        install=GH + "xiaoyulu95/SSMD@b90d2541d8c5f76733c4ea151972f237c126ffe4;cran::nnls",
        source="GitHub:xiaoyulu95/SSMD", version="0.1.0", pin="b90d2541d8c5f76733c4ea151972f237c126ffe4",
        license="none (placeholder 'What license is it under?')", lic_note="DESCRIPTION template left unfilled",
        expr_ev="not stated; ran with TPM", registration="not required",
        notes="use the ORIGINAL repository (the zy26/SSMD fork lacks R2_two_mat); needs the undeclared package nnls; 6 samples or more (3 and 2 fail); 90-112 s per run; rows are data-dependent marker modules (21-30), not fixed cell types"),
    "immucc": dict(species="mouse", status="excluded", min_n="", single="", indep="", cohort="",
        needs="", options="", cls="fraction", packages="", install="",
        source="GitHub:wuaipinglab/ImmuCC", version="", pin="", license="not stated", lic_note="README",
        expr_ev="", registration="required",
        notes="the algorithm is CIBERSORT's SVR and the script tells users to request CIBERSORT from cibersort.stanford.edu; its tissue signature matrices (CSV) are open and can be reused with an open solver"),

    "seqimmucc": dict(species="mouse", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="fraction", packages="", install="", engine="native",
        source="native port of seq-ImmuCC LLSR (Chen et al. 2018); signature matrix from immunedeconv@" + DATA_SHA[:7],
        version="bulkdeconv", pin="data:" + DATA_SHA, license="GPL-3 (bulkdeconv)",
        lic_note="port of code marked as provided by the seq-ImmuCC authors; the signature matrix has no stated licence (inst/COPYRIGHTS)",
        expr_ev="raw counts (paper and immunedeconv documentation); the code log2-transforms values above 50; a TPM-like matrix also ran, but only counts are accepted here",
        registration="not required", expr="count", n_out="10", secs="0.1",
        notes="LLSR only (the SVR variant needs CIBERSORT and is not offered); 10 immune cell types in percent summing to 100; the original computes a quantile normalisation it never uses, so none is applied; values identical alone, in n=5 and in n=20"),
    "base": dict(species="mouse", status="conditional", min_n="2", single="unsupported", indep="no (diff 0.65 between n=6 and n=3)", cohort="yes",
        needs="n_permutations (default 100), seed (default 1)", options="n_permutations|seed", cls="enrichment_score", packages="", install="", engine="native",
        source="native port of BASE (Varn et al. 2017) with the mouse compendium from immunedeconv@" + DATA_SHA[:7],
        version="bulkdeconv", pin="data:" + DATA_SHA, license="GPL-3 (bulkdeconv)",
        lic_note="port of code adapted from Varn et al.; the compendium has no stated licence and is stored rounded to 6 decimals (inst/COPYRIGHTS)",
        expr_ev="TPM (immunedeconv documentation)", registration="not required", n_out="15", secs="51",
        notes="51 s for 6 samples and 57 s for 20 samples at 100 permutations (the original code: 175 s for 6); quantile normalisation and per-gene median centring make scores relative to the cohort; scores are z-scored over samples before sub-types are combined, so 1 sample is not possible; permutations are seeded (option seed); 15 mouse cell types"),
    "abis": dict(species="human", status="ready", min_n="1", single="native", indep="yes (diff 0)", cohort="no",
        needs="", options="", cls="abundance_estimate", packages="MASS", install="", engine="native",
        source="native port of ABIS (Monaco et al. 2019), RNA-seq signature from immunedeconv@" + DATA_SHA[:7],
        version="bulkdeconv", pin="data:" + DATA_SHA, license="GPL-3 (bulkdeconv)",
        lic_note="port of the ABIS Shiny application code; the signature matrix has no stated licence (inst/COPYRIGHTS)",
        expr_ev="TPM (immunedeconv documentation; RNA-seq signature)", registration="not required", n_out="17", secs="0.1",
        notes="signature built from sorted blood (PBMC) cell types, not from tumour tissue: 17 blood cell subsets as robust-regression coefficients times 100; uses MASS (a recommended package that most R installations include; installed into the private library when missing)"),
    "timer": dict(species="human", status="conditional", min_n="2", single="unsupported", indep="no (diff 0.04 between n=4 and n=3)", cohort="yes",
        needs="cancer_type (32 TCGA types, required)", options="cancer_type", cls="abundance_score", packages="sva", install="", engine="native",
        source="native port of TIMER (Li et al. 2016); reference data from immunedeconv@" + DATA_SHA[:7],
        version="bulkdeconv", pin="data:" + DATA_SHA, license="GPL-3 (bulkdeconv)",
        lic_note="port of TIMER code (GPL-2.0); the curated immune reference derives from public sorted-cell expression data (inst/COPYRIGHTS)",
        expr_ev="TPM (immunedeconv documentation)", registration="not required", n_out="6", secs="5.7",
        notes="signatures are cancer-type specific (a different type changed the values by up to 0.02); batch correction against a reference (ComBat, package sva) makes values depend on the other samples; needs at least 2 samples; 6 cell types (B, CD4 T, CD8 T, neutrophil, macrophage, DC)"),
}
order = ["quantiseq", "mcp_counter", "xcell", "xcell2", "estimate", "consensus_tme", "immucellai2", "imsig",
         "spatialdecon_safetme", "epic", "immucellai", "abis", "timer", "mmcp_counter", "dcq", "immucellai_mouse", "ssmd", "seqimmucc", "base", "immucc"]


def measured(method):
    if "n_out" in M[method]:
        return M[method]["secs"], M[method]["n_out"]
    key = "xcell2[TMECompendium]" if method == "xcell2" else method
    r = multi.get(key)
    if method == "ssmd":
        return "90-112", "21-30"
    if not r or r["status"] != "ok":
        return "", ""
    return r["seconds"], r["out_dim"].split("x")[-1]


fields = ["method", "species", "status", "expr_accepted", "min_samples_n", "single_sample", "sample_independent",
          "cohort_dependent", "engine", "needs", "options", "output_class", "comparability", "n_outputs_measured", "runtime_s",
          "packages", "install", "source", "version", "pin", "license", "license_note", "registration",
          "expr_evidence", "validated_on", "notes"]
with open(out_dir / "methods.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    for m in order:
        c = M[m]
        secs, n_out = measured(m)
        usable = c["status"] in ("ready", "conditional")
        w.writerow(dict(method=m, species=c["species"], status=c["status"], expr_accepted=c.get("expr", "tpm") if usable else "",
                        min_samples_n=c["min_n"], single_sample=c["single"], sample_independent=c["indep"],
                        cohort_dependent=c["cohort"], engine=c.get("engine", "external"), needs=c["needs"], options=c["options"], output_class=c["cls"],
                        comparability=CLASS_NOTE[c["cls"]], n_outputs_measured=n_out, runtime_s=secs,
                        packages=c["packages"], install=c["install"], source=c["source"], version=c["version"],
                        pin=c["pin"], license=c["license"], license_note=c["lic_note"],
                        registration=c["registration"], expr_evidence=c["expr_ev"],
                        validated_on=VALIDATED if usable else "", notes=c["notes"]))

shared = {"TabulaSapiensBlood": "0.75", "MouseRNAseqData": "0.83", "ImmGenData": "0.85", "TabulaMurisBlood": "0.88"}
presets = {  # name: (species, platform, tissue, default)
    "TMECompendium": ("human", "RNA-seq (TPM)", "tumour", "yes"),
    "PanCancer": ("human", "scRNA-seq (counts)", "tumour", ""),
    "ImmuneCompendium": ("human", "RNA-seq (TPM)", "immune/blood", ""),
    "BlueprintEncode": ("human", "RNA-seq (TPM)", "mixed", ""),
    "LM22": ("human", "microarray (RMA)", "mixed; LM22 originates from CIBERSORT (non-commercial)", ""),
    "TabulaSapiensBlood": ("human", "scRNA-seq (counts)", "blood, lymph node, spleen, thymus, bone marrow", ""),
    "MouseRNAseqData": ("mouse", "RNA-seq (TPM)", "mixed", "yes"),
    "ImmGenData": ("mouse", "microarray (RMA)", "immune/blood", ""),
    "TabulaMurisBlood": ("mouse", "scRNA-seq (counts)", "bone marrow, spleen, thymus", ""),
}
with open(out_dir / "presets.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["method", "preset", "species", "platform", "tissue", "default", "n_outputs_measured", "runtime_s",
                "shared_gene_fraction_on_test_data", "run_status_multi", "run_status_duplicated"])
    for name, (sp, plat, tis, dflt) in presets.items():
        key = f"xcell2[{name}]"
        m, d = multi.get(key, {}), dup.get(key, {})
        w.writerow(["xcell2", name, sp, plat, tis, dflt, m.get("out_dim", "").split("x")[-1], m.get("seconds", ""),
                    shared.get(name, ">=0.9"), m.get("status", ""), d.get("status", "")])
print("wrote", out_dir / "methods.csv", "and", out_dir / "presets.csv")
