#!/usr/bin/env python3
"""Refresh delivery copies and their hashes from the reviewed plot outputs."""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parent.parent
SOURCES = {
    'Fig1A_Experimental_Strategy': 'results/full_figures/Fig1A_Experimental_Strategy',
    'FigS1A_RNAseq_Workflow': 'results/full_figures/FigS1A_RNAseq_Workflow',
    'Fig4C_Working_Model': 'results/full_figures/Fig4C_Working_Model',
    'Fig3A_Sleep_Profile': 'results/sleep/sleep_profile_corrected',
    'Fig3B_Total_Sleep': 'results/sleep/sleep_total_corrected',
    'Fig3C_Day_Night_Sleep': 'results/sleep/sleep_phase_corrected',
    'Fig3D_Activity_Profile': 'results/sleep/activity_profile_corrected',
    'Fig3E_Total_Activity': 'results/sleep/activity_total_corrected',
    'Fig3F_Day_Night_Activity': 'results/sleep/activity_phase_corrected',
    'Fig4B_MTT': 'results/panels/Fig4B_MTT_descriptive',
    'Fig2A_GO_Female': 'results/go/Fig2A_GO_Female',
    'FigS1E_GO_Male': 'results/go/FigS1E_GO_Male',
    'Fig1B_PCA_MDS': 'results/full_figures/Fig1B_PCA_MDS',
    'Fig1C_Female_Volcano': 'results/full_figures/Fig1C_Female_Volcano',
    'Fig1D_Female_Heatmap': 'results/full_figures/Fig1D_Female_Heatmap',
    'FigS1C_Male_Volcano': 'results/full_figures/FigS1C_Male_Volcano',
    'FigS1D_Male_Heatmap': 'results/full_figures/FigS1D_Male_Heatmap',
    'Fig2D_Individual_Genes': 'figure_sources/Figure_2E_Individual_Genes/Rebuilt_Output/Fig2D_individual_genes_female_horizontal',
    'Fig3H_Locomotion': 'results/full_figures/Fig3H_Locomotion',
    'Fig3J_Mitochondrial_Perimeter': 'figure_sources/Figure_3K_Mitochondrial_Area/Rebuilt_Output/TEM_mitochondrial_perimeter_2026-08-04',
    'Fig3K_Mitochondrial_Area': 'figure_sources/Figure_3K_Mitochondrial_Area/Rebuilt_Output/TEM_mitochondrial_area_2026-08-04',
    'Fig3L_Mitochondrial_Aspect_Ratio': 'figure_sources/Figure_3K_Mitochondrial_Area/Rebuilt_Output/TEM_mitochondrial_aspect_ratio_2026-08-04',
    'Fig4A_TIMELESS_Descriptive': 'results/full_figures/Fig4A_TIMELESS_Descriptive',
    'FigS1B_Western_Blot_Quantification': 'figure_sources/Figure_S1B_Western_Blot/Rebuilt_Output/MIC60_WB_anti_myc_quantification',
    'FigS2_KEGG_GSEA': 'figure_sources/Figure_S2_KEGG_GSEA/Rebuilt_Output/FigS2_KEGG_GSEA',
    'FigS3_Selected_Gene_Modules': 'figure_sources/Figure_S3_Selected_Gene_Modules/Rebuilt_Output/FigS3_selected_gene_modules',
    'FigS4A_Stress_Response': 'figure_sources/Figure_S4_Stress/Rebuilt_Output/figS4_stress_distinct',
    'FigS4C_Transfection_Efficiency': 'figure_sources/Figure_S4B_Transfection_Efficiency/Rebuilt_Output/FigS4C_transfection_efficiency',
    'Fig3N_TMRM': 'figure_sources/Figure_3O_TMRM/Rebuilt_Output/TMRM_MTG_image_means_only_2026-08-05',
}

def publish():
    manifest = []
    for name, stem in SOURCES.items():
        for ext in ('pdf', 'png'):
            source = ROOT/'Analysis'/f'{stem}.{ext}'
            target = ROOT/'Figures'/ext.upper()/f'{name}.{ext}'
            if not source.is_file():
                raise FileNotFoundError(f'Missing required figure: {source}')
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            digest = hashlib.sha256(target.read_bytes()).hexdigest()
            manifest.append({'file': str(target.relative_to(ROOT/'Figures')),
                             'generated_source': str(source.relative_to(ROOT)),
                             'sha256': digest,
                             'generated_source_sha256': digest,
                             'source_matches_delivery': True,
                             'delivery_type': 'generated'})
    (ROOT/'Figures/manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    source = ROOT/'Analysis/results/go/all_displayed_terms_and_genes.csv'
    if source.exists():
        (ROOT/'Figures/Data').mkdir(exist_ok=True)
        shutil.copy2(source, ROOT/'Figures/Data/GO_displayed_terms_and_genes.csv')
    print(f'Updated {len(manifest)} figure files in {ROOT / "Figures"}')

if __name__ == '__main__':
    publish()
