# scsa
Sperm Chromatin Structure Analysis (SCSA) program for R

## Overview
SCSA is an assay used to test for DNA damage in sperm cells using flow cytometry. Details of the protocol are given here:

Evenson D, Jost L, Gandour D et al. (1995) Comparative sperm chromatin structure assay measurements on epiillumination and orthogonal axes flow cytometers. Cytometry, 19, 295–303.

## Installation

```
install.packages("scsa")
```

## Usage
This package takes raw data from the flow cytometer and calculates the SCSA statistics used to assess DNA damage in the sample. The input is a matrix with the "green" channel in the rows and the "red" channel in columns. Each cell in the matrix contains the count of biological cells that emitted that combination of wavelengths in the "green" and "red" channel while passing through the flow cytometer. The matrix must have the same number of rows and columns, one for each channel. 

Consult your flow cytometer software manual for how to export the raw data. 

To analyze your flow cytometry data:

1) Import the raw data from a sample using read.csv()

2) run scsa() on the imported data.

The program will run you through a series of histograms to select the values for Alphat and HDS. Calling plot() generates a figure:

![SCSA_plot](https://github.com/chadsmith123/scsa/blob/master/exdata/sample_plot.png)

Thats it! The program will guide you through the analysis steps and return an scsa object with the statistics.
