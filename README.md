# Organic Photovoltaic Predictor 
This calculator is published at (http://info.eecs.northwestern.edu/OPVPredictor)

## Pre-requisites: 
1. Python 2.7 (or higher) 
2. Sklearn 0.14 (or higher) 
3. Rdkit 2012.9 (or higher)  
4. Numpy 1.4.1 (or higher) 

## List of Files: 
* \__init__.py: Core python file that launches the predictor  
* atomModel.pkl: Pickle file containing the model for predicting HOMO from Atom Pair fingerprints 
* maccsModel.pkl: Pickle file containing the model for predicting HOMO from MACCS fingerprints 
* static directory: This directory contains the static elements (boostrap theme files, images etc)
* templates directory: This directory contains 2 files 
  - index.html: HTML file for the homepage (request page)
  - OPV.html: HTML file for the response page 


## Developer Team

The code was developed by the <a href="http://cucis.ece.northwestern.edu/">CUCIS</a> group at the Electrical and Computer Engineering Department at Northwestern University. 

1. Arindam Paul (arindam.paul@eecs.northwestern.edu)
2. Ankit Agrawal (ankitag@eecs.northwestern.edu)
3. Wei-keng Liao (wkliao@eecs.northwestern.edu)
4. Alok Choudhary (choudhar@eecs.northwestern.edu)


The development team would like thank our collaborator <a href="https://www.feinberg.northwestern.edu/faculty-profiles/az/profile.html?xid=40386">Prof. Alona Furmanchuk</a> from <a href="https://www.feinberg.northwestern.edu/">Northwestern Feinberg School of Medicine</a>. 


## Citation
If you use this code or data, please cite:

Arindam Paul, Alona Furmanchuk, Wei-keng Liao, Alok Choudhary, Ankit Agrawal. Property Prediction of Organic Donor Molecules for Photovoltaic Applications using Extremely Randomized Trees. Journal of Molecular Informatics, 2019


## Questions/Comments:

email: arindam.paul@eecs.northwestern.edu or ankitag@eecs.northwestern.edu</br>
Copyright (C) 2019, Northwestern University.<br/>
See COPYRIGHT notice in top-level directory.

## Funding Support

This work was performed under the following financial assistance awards 70NANB14H012 and 70NANB19H005 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD). Partial support is also acknowledged from DOE awards DE-SC0014330, DE-SC0019358.
