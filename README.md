# Encore 

Encore allows scientists to test hypotheses using large scale sequencing data
while hiding many of the complicated technical issues associated with working
on terabyte scale data.  In addition to performing the most common association
tests for single variants and groups of rare variants, Encore enables users to
annotate, share and interact with results.  Access to Encore is controlled by a
white list of users who can upload phenotypes and request analyses.  All
per-variant and per-gene summary statistics are returned.  To maintain
confidentiality of research participant data, individual level data cannot be
downloaded.  By optimizing the most common analyses and providing rich ways to
interact with results, Encore offers an exciting platform to investigate
genetic hypotheses.

# Getting Started

To get access to Encore, please contact Snehal Patil at snehal@umich.edu to be added to the whitelist. You will also need to log into the UM VPN while using Encore. 

## Preparing Your Phenotype File

Encore already contains the genotypes on the backend, so you do not need to upload or directly work with any genetic data files. However, you will need to upload a file containing the phenotypes you wish to analyze. This file should also contain any covariates that you wish to adjust for in your model.

Click "Upload Phenotype File" to upload a file for analysis.

Please ensure your phenotype file meets the following specifications:
* CSV or tab-delimited plain text
* Each row contains same number of fields
* Column headers included as first row
* Contains at least two columns: sample ID and phenotype values

Example Phenotype File:
```
ID                       PHENOTYPE      AGE       CHIP
XXXX-XXXX-XXXX-0001           0         45        1
XXXX-XXXX-XXXX-0002           1         32        1
XXXX-XXXX-XXXX-0003           0         76        0
XXXX-XXXX-XXXX-0004           0         55        1
```

IMPORTANT: Sample ID MUST be `DeID ID` from DataDirect. MRNs or biobank IDs will not work with Encore.

## Running Your GWAS

Once your phenotype file has been successfully uploaded, you can start submitting GWAS analysis jobs. Once submitted, jobs will be scheduled on the Great Lakes computing cluster and run in the order they were received as compute resources become available. You can check the status of your job at any time under the "Jobs" tab in the "Status" column.

To prepare and submit a GWAS job, click on "Build Model" on the pop-up that appears directly after uploading your phenotype file. Or, you can click on "Create a New Job from Existing Phenotypes" on the main page.

`Genotypes`: By default, "MGI_Freeze3" is selected. This dataset is TOPMed-imputed and on build 38. Currently this is the only set of genotype data available for analysis. As we continue to add new freezes, you will find more options here.

`Phenotypes`: Select a phenotype file that you have uploaded.

`Response`: Select the column name for your main response variable/outcome.

`Covariates`: Select any covariates that you wish to adjust for in your model, such as age, sex, etc.

`Pop adjust`: Select genetic principal components to adjust for population stratification in your model. We recommend using at least 5, but modern GWAS analyses typically use up to 10 PCs. These PCs are calculated in relation to the entire cohort. If you wish to use a different set of PCs (e.g. computed only in Europeans or only in your analysis sample), you will need to upload them as part of your phenotype file and select them under the `Covariates` input field.

`Model`: Encore currently supports a simple Linear Wald test for quantitative traits (EPACTS) and mixed model analyses for both quantitative and binary traits (SAIGE). We recommend SAIGE as it adjusts for relatedness among samples as well as case-control imbalance for binary traits. If you have not removed related individuals from your dataset, please do not use the EPACTS implementation. Currently these are the only types of models available in Encore, but we will make more available in the future. If you are interested in running a test that is currently not supported, please contact us.


# Developer Guide

## Installing python modules

To install all the required python modules into a virtal environment

     python3 -m venv venv/
     source venv/bin/activate
     pip install -r requirements.txt

## Run dev instance

To run a simple dev instance of Encore, you can run

    ./rundev.sh

## Apache Configuration

Currently Encore is deployed using WSGI with apache. You can 
install apache and the WSGI module with

    apt-get install apache2 python-setuptools libapache2-mod-wsgi-py3 python3-dev libmysqlclient-dev

A sample wsgi file is included at encore.wsgi.example. You should
copy this file to encore.wsgi and fill in the correct path
for your server.

A sample apache conf file is included at encore.conf.example. You should
copy this file to /etc/apache2/sites-available/encore.conf and
fill in the correct URL and path for your server.

You can enable the configuration with

    sudo a2ensite encore
    sudo systemctl restart apache2

## Email configuration

Encore will try to send email via localhost. You can set up
postfix to send or redirect these emails how you like. You can
install postfix with

    sudo apt-get install postfix

A simple installation would just choose "Internet Site" and set

    inet_interfaces = localhost

in the `/etc/postfix/main.cf` config file.



## Building Executable tools

     mkdir build
     cd build
     cmake ..
     make
