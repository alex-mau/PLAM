Aug, 2011

This subdirectory contains the following Matlab source codes:

selector.m            The main code that solves for Dantzig selectors.

explicitMatrix.m      This code creates function handles for selector.m.

spg2.m                This code solves the ADM subproblem by spectral projected gradient method.

linesearch.m          This code does the Armijo linesearch for spg2.m.

func.m                This code computes the function value for spg2.m.

grad.m                This code computes the gradient for spg2.m.

soft_thresh.m         This code computes the soft thresholding of L-1 norm for spg2.m.

postprocess.m         This code computes the two stage Dantzig selector from the output of selector.m.

Implementation and numerical experience with selector.m are described in the paper: 
    Zhaosong Lu, Ting Kei Pong and Yong Zhang
    "An Alternating Direction Methods for Finding Dantzig Selectors",
    Submitted.
This code was last updated on August 21, 2011.

Questions/comments/suggestions about the codes are welcome.  

Zhaosong Lu,   zhaosong@sfu.ca
Ting Kei Pong, tkpong@math.washington.edu
Yong Zhang,    yza30@sfu.ca