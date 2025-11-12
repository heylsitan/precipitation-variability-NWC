# Precipitation Variability Analysis for Northwest China

**Associated Paper**: [To be updated with your JGRA manuscript title and DOI]
**Authors**: [Your names]
**Journal**: Journal of Geophysical Research: Atmospheres (JGRA)
**Contact**: heyongli@lzu.edu.cn

This repository contains the analysis code for studying precipitation variability patterns in Northwest China, including internal vs external variability, trend analysis, and moisture source tracking.

## Project Structure

The code has been reorganized from the original Chinese version with English comments and improved file naming while preserving all original functionality and directory structure.

```
organized_project/
├── station_merra2_comparison/           # Station-MERRA2 Comparison Analysis
├── variability_analysis/                # Internal/External Variability Analysis
├── spatial_climatology/                 # Spatial Climatology and Mapping
├── advanced_analysis/                   # Advanced Analysis (SPI, Moisture Tracking, Drought)
├── config/                              # Configuration files
├── docs/                                # Documentation
├── data/                                # Input data directory
└── outputs/                             # Output figures directory
```

## File Naming Convention

Files maintain the original `figureX_stepX` format but with meaningful English descriptions:
- Format: `figure[number]_[descriptive_name].py` for main scripts
- Format: `figure[number]_step[step_number]_[descriptive_name].py` for step scripts
- Format: `figure[number]_step[step_number]_[descriptive_name].sh` for shell scripts

## Analysis Components by Directory

### Station-MERRA2 Comparison Analysis
- **figure1_station_merra2_comparison.py** - Main comparison between station and MERRA-2 data
- **figure1_step1_station_data_loading.py** - Station data loading and processing
- **figure1_step2_merra2_data_processing.py** - MERRA-2 reanalysis data processing
- **figure1_step3_data_download.sh** - Data download and merging script
- **figure1_step4_trend_analysis.py** - Trend analysis for precipitation data
- **figure1_step5_correlation_computation.py** - Correlation computation
- **figure1_step6_visualization_setup.py** - Visualization setup with spatial analysis

### Internal/External Variability Analysis
- **figure3_internal_external_variability_analysis.py** - Main variability decomposition analysis
- **figure3_step1_variability_data_loading.py** - Variability data loading step
- **figure3_step2_trend_computation.py** - Trend computation for variability components
- **figure3_step3_spatial_data_processing.py** - Spatial data processing step
- **figure3_step4_mapping_visualization.py** - Mapping and visualization step

### Spatial Climatology and Mapping
- **figure4_moisture_source_regional_mapping.py** - Moisture source regional mapping
- **figure5_precipitation_climatology_mapping.py** - Precipitation climatology spatial mapping
- **figure6_trend_significance_mapping.py** - Trend significance spatial mapping
- **figure7_spatial_trend_analysis.py** - Spatial trend analysis
- **figure7_step1_climatology_data_processing.py** - Climatology data processing step

### Advanced Analysis
- **figure8_standardized_precipitation_index.py** - SPI calculation and drought analysis
- **figure9_moisture_tracking_analysis.py** - Comprehensive moisture source tracking
- **figure10_trend_comparison_analysis.py** - Atmospheric circulation trend comparison
- **figure11_spatial_spi_analysis.py** - Spatial Standardized Precipitation Index analysis
- **figure12_regional_precipitation_analysis.py** - Regional precipitation and circulation patterns
- **figure9_step1_moisture_data_preprocessing.py** - Moisture data preprocessing step
- **figure11_step1_spi_data_preprocessing.py** - SPI data preprocessing step

## Requirements

- Python 3.7+
- Required packages: numpy, pandas, matplotlib, xarray, cartopy, scipy, netCDF4, cmaps, geopandas, scikit-learn, shapely
- Additional tools: CDO (Climate Data Operators) for some shell scripts

## Usage

1. Update the configuration file: `config/config.py`
2. Place your data in the appropriate directories
3. Run individual analysis scripts as needed
4. Follow the step-by-step workflow for complex analyses

## Data Requirements

- Station precipitation data (Excel format)
- MERRA-2 reanalysis data (NetCDF format)
- Shapefiles for regional boundaries
- Moisture tracking results (NetCDF format)

## Key Features Maintained

- **Original Code Logic**: All scientific calculations and algorithms preserved exactly
- **Step-by-Step Workflow**: Maintains the original research methodology
- **Chinese to English Translation**: All comments translated for international accessibility
- **Professional Naming**: Descriptive English names while preserving figure numbering
- **Academic Standards**: Suitable for publication and collaboration

## Geographic Focus

- **Region**: Northwest China (70-115°E, 30-51°N)
- **Time Period**: 1982-2020 (39 years)
- **Data Sources**: Station observations, MERRA-2 reanalysis

## Analysis Capabilities

- **Precipitation trend analysis** with statistical significance testing
- **Internal vs external variability decomposition**
- **Moisture source tracking and regional contribution analysis**
- **Standardized Precipitation Index (SPI) calculation for drought identification**
- **Spatial climatology and trend mapping**
- **Atmospheric circulation pattern analysis**

## License

This code is provided for academic research purposes under MIT License.

## Note

This is an English-translated version of the original Chinese analysis code, with improved file naming and documentation while preserving all original functionality and directory structure. The code maintains the original scientific methodology and is ready for international collaboration and publication.