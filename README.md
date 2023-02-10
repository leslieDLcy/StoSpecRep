# StoSpecRep

As the name suggests, I'm reorganizing the `spectral representation` package. 

## Usage
It mainly contains several parts:
1. Formulating PSD models;
2. Estimating EPSD or PSD from realizations;

## Note
By default, given a recording, we assume it as the target and hence compute the EPSD of this recording. But also, the `CWTx` class can be used for an ensemble of recordings as the *externaldata*, where we are computing an ensemble of EPSDs. 