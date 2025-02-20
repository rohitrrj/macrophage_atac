# Macrophage ATAC-seq Visualization App
![Docker](https://img.shields.io/badge/docker-ready-blue.svg)
![R-Shiny](https://img.shields.io/badge/R--Shiny-app-brightgreen.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

An interactive R Shiny application for visualizing and analyzing ATAC-seq data from macrophage samples. This application provides a user-friendly interface to explore chromatin accessibility patterns in macrophages.

## Table of Contents
- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Data Description](#data-description)
- [Docker Deployment](#docker-deployment)
- [Configuration](#configuration)
- [Contributing](#contributing)

## Features
- Interactive visualization of ATAC-seq peaks
- Analysis of differential accessibility between conditions
- Integration with gene expression data
- Dynamic filtering and sorting capabilities
- Customizable visualization parameters
- Export functionality for plots and data

## Quick Start

Using Docker:
```bash
docker pull yourusername/macrophage-atac
docker run -p 3838:3838 yourusername/macrophage-atac
```

The app will be available at: http://localhost:3838

## Installation

### Local Installation
1. Clone the repository:
```bash
git clone https://github.com/yourusername/macrophage_atac.git
cd macrophage_atac
```

2. Install R dependencies:
```R
install.packages(c("shiny", "ggplot2", "DT", "plotly"))
```

3. Run the app:
```R
shiny::runApp("app")
```

### Docker Installation
Build the Docker image:
```bash
docker build -t macrophage-atac .
```

## Usage

1. Launch the application
2. Upload or select pre-loaded data
3. Use the sidebar controls to:
   - Select samples/conditions
   - Adjust visualization parameters
   - Filter data points
4. Interact with plots and tables
5. Export results as needed

## Data Description

The application works with:
- Processed ATAC-seq peak files
- VST-normalized count matrix (`vstNormalizedCounts_Macrophage.txt`)
- Gene annotations (`Gene_Symbols.txt`)
- Differential accessibility results (`HCvsCAD.txt`)

### Data Format Requirements

Peak Files:
```
chromosome start end signal
chr1    1000    1500    5.2
chr1    2000    2500    3.1
...
```

Count Matrix:
```
GeneID Sample1 Sample2 Sample3
GENE1  10.5    11.2    9.8
GENE2  8.7     7.9     8.1
...
```

## Docker Deployment

The application is containerized using Docker for easy deployment:

```dockerfile
FROM rocker/shiny:4.1.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev

# Copy application files
COPY app /srv/shiny-server/app
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

# Set permissions
RUN chmod +x /usr/bin/shiny-server.sh

# Expose port
EXPOSE 3838

# Start Shiny server
CMD ["/usr/bin/shiny-server.sh"]
```

## Configuration

### Shiny Server Configuration
The application uses a custom Shiny server configuration (`shiny-customized.config`):
```
run_as shiny;
preserve_logs true;
access_log /var/log/shiny-server/access.log tiny;
server {
  listen 3838;
  location / {
    site_dir /srv/shiny-server/app;
    log_dir /var/log/shiny-server;
    directory_index on;
  }
}
```

### Application Settings
Modify `app/app.R` for:
- Maximum upload size
- Cache settings
- Plot parameters
- Default filters

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

### Development Guidelines
- Follow R code style guidelines
- Add tests for new features
- Update documentation
- Test Docker deployment

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this application in your research, please cite:
```
Author et al. (Year). Macrophage ATAC-seq Visualization: An interactive tool for analyzing chromatin accessibility.
```

## Acknowledgments
- Data processing pipeline: [ATACseq_Pipeline](../ATACseq_Pipeline)
- R Shiny framework
- Supporting institutions and funding
