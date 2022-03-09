## Constant-Q transform toolbox for music processing (2010)

This is the reference implementation accompanying the paper

Schörkhuber, Christian, and Klapuri, Anssi. "[Constant-Q transform toolbox for music processing](https://iem.kug.ac.at/fileadmin/media/iem/projects/2010/smc10_schoerkhuber.pdf)." In 7th sound and music computing conference, Barcelona, Spain, pp. 3-64. 2010.

**NOTE**: Since the site that hosted the download link to this toolbox is offline, this repository contains the original Matlab code accompying the above pubication. **The code is not actively maintained.**

### Abstract:
This paper proposes a computationally efficient method for computing the constant-Q transform (CQT) of a time-domain signal. CQT refers to a time-frequency representation where the frequency bins are geometrically spaced and the Q-factors (ratios of the center frequencies to bandwidths) of all bins are equal. An inverse transform is proposed which enables a reasonable-quality (around 55dB signal-to-noise ratio) reconstruction of the original signal from its CQT coefficients. Here CQTs with high Q-factors, equivalent to 12--96 bins per octave, are of particular interest. The proposed method is flexible with regard to the number of bins per octave, the applied window function, and the Q-factor, and is particularly suitable for the analysis of music signals. A reference implementation of the proposed methods is published as a Matlab toolbox. The toolbox includes user-interface tools that facilitate spectral data visualization and the indexing and working with the data structure produced by the CQT.
