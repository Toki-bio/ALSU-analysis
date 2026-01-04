
https://toki-bio.github.io/ALSU-analysis/steps/step1.html
# ALSU Genotyping Pipeline - Workflow Documentation

## Project Overview

This directory contains comprehensive HTML documentation for the ALSU Uzbek Pregnancy Loss Cohort genotyping imputation and quality control pipeline.

## Structure

- **index.html** - Main page with roadmap of all 6 pipeline steps
- **steps/** - Individual step documentation pages (step1.html through step6.html)
- **daily/** - Daily activity logs with terminal output and timestamps
- **config.json** - Configuration file for workserver connection
- **SETUP_GUIDE.html** - Complete setup and maintenance guide

## Quick Start

1. Open `index.html` in a web browser to view the main documentation
2. Click on any step card to view detailed documentation
3. Use the "Daily Logs" tab to view daily activity and terminal output
4. See `SETUP_GUIDE.html` for instructions on setting up workserver integration

## Key Features

### Main Page (index.html)
- Roadmap view with all 6 pipeline steps
- Quick metrics showing:
  - 1,074 final samples
  - 10.8 million variants
  - 6 sequential processing steps
  - December 2025 completion date

### Step Pages (steps/step1.html - step6.html)
Each step has 3 tabs:
1. **Overview** - Purpose, key metrics, rationale
2. **Technical Details** - Parameters, commands, files, quality checks
3. **Chronological Log** - Timeline of events with timestamps

Step breakdown:
1. Sample Missingness Filter (Dec 15)
2. IBD Deduplication (Dec 16)
3. SNP QC & VCF Export (Dec 17)
4. Michigan Imputation Server (Dec 18-22)
5. Sample ID Normalization - Homoglyph Fix (Dec 23)
6. Post-Imputation QC & Finalization (Dec 26)

### Daily Logs (daily/YYYY-MM-DD.html)
Terminal-style documentation with:
- Session summary box
- Timestamped log entries
- Complete command output
- Status indicators

## Configuration

Edit `config.json` to configure:
- Workserver hostname and credentials
- SSH connection settings
- Base directory on remote server
- Log sync interval

## Workserver Integration

To set up automatic log syncing from Biotech2024:

1. Update credentials in config.json
2. Set up SSH key authentication (recommended)
3. Create sync script using provided template in SETUP_GUIDE.html
4. Schedule with Windows Task Scheduler or cron job

## File Naming Convention

- Daily logs: `daily/YYYY-MM-DD.html` (e.g., `daily/2025-12-15.html`)
- Multi-day logs: `daily/YYYY-MM-DD-DD.html` (e.g., `daily/2025-12-18-22.html`)
- Step pages: `steps/step#.html`

## Original Documentation

Original markdown documentation archived in `old_logs/`:
- complete_pipeline_summary.md
- alsu_genotype_imputation_pipeline_full_technical_log.md
- step6_technical_log.md
- ConvSK_QC_and_VCF_report.md

## Maintenance

- Update "Last updated" date in footer when modifying pages
- Review links quarterly for broken references
- Archive old documentation as pipeline progresses
- Keep config.json credentials secure

## Support & Updates

For setup assistance or questions about the documentation system, refer to:
- SETUP_GUIDE.html - Complete setup instructions
- Existing step pages - Examples of documentation structure
- Daily log templates - Format for new activity logs

---

**Project**: ALSU Uzbek Pregnancy Loss Cohort  
**Status**: ✓ Completed (Ready for Population Genetics Analyses)  
**Final Dataset**: UZB_imputed_HQ_clean (1,074 samples, 10.8M variants)  
**Last Updated**: January 2, 2026
