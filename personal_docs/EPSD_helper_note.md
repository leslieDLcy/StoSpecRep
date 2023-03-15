Here below we have some helpful links to recap:

1.[Practical introdution to CWT](https://uk.mathworks.com/help/wavelet/ug/practical-introduction-to-continuous-analysis.html#PracIntroCWAnalysisExample-7)
2.[Continuous 1D CWT matlab](https://uk.mathworks.com/help/wavelet/ref/cwt.html?s_tid=srchtitle#bvb0t8f-1_1)
3. [Wavelet Transform in Matlab](https://uk.mathworks.com/discovery/wavelet-transforms.html)


## Personal Note
By default, given a recording, we assume it as the target and hence compute the EPSD of this recording. But also, the `CWTx` class can be used for an ensemble of recordings as the *externaldata*, where we are computing an ensemble of EPSDs. 


```
Note:
Months ago, we first defined a 'SRM' class which turns out can only be used with 
spectrogram obtained by STFT, but not with scalogram obtained by Wavelet Transform.
The reason is the shape of EPSD estimates, STFT give very small shape, eg. (57, 129)
which needs to be interpolated in spectrum to further generate sample realizations.

But Wavelet Transform can give a large Swt shape such that we don't need to interpolate anymore.
Therefore, it just needs a `SRM_formulat` function (see below).
```