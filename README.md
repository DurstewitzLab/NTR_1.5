#Neural Trajectories Reconstruction Matlab Toolbox 1.5 beta

# I. Overview
This toolbox performs statistical analyses on nonlinear time-series of multivariate
neural responses –of different kinds- during a cognitive task. 
Therefore it is assumed that the separate cognitive epoch's tasks are observable (e.g. a stimulus presentation, a 
successful choice, reward acquisition, movement to a definite place like an arm) and are
labelled accordingly. This requirement will be relaxed in future versions.

It is also assumed that all multivariate neural responses are time-ordered and
simultaneously recorded. Convergence or divergence of neural responses trajectories will
be evaluated in those states defined by the distinct time-windows corresponding to
cognitively relevant epochs of the task. Statistical analyses will indicate the probability of
those states to “behave” like attracting regions of neural responses (or rather they will
suggest more likely transient dynamics). This toolbox is intended to be used in different
types of neurophysiological response. Nevertheless it was developed and testes only on
electrophysiological data (not on BOLD or EEG responses at least in the present
moment).

An optimal reconstruction of neural trajectories properties may require high-dim,
nonlinear neural interactions spaces. For that reason, kernel methods will be used in the
statistical analyses of such state spaces.

A demonstrative dataset is provided, please type “>ntr;”. Besides in Balaguer et
al. (2011), a brief overview of such techniques can be found in [Durstewitz and Balaguer
(2010), Statistical Approaches for Reconstructing Neuro-Cognitive Dynamics from High-
Dimensional Neural Recordings. Neuroforum 4/10: 266-276](https://www.degruyter.com/view/j/nf.2010.16.issue-4/s13295-010-0011-0/s13295-010-0011-0.xml?format=INT)
 and in references therein (e.g. Chuckland et al., 1999 etc.). Automatic “data-smoothing” methods (e.g. for
estimate firing rates from spike trains) are not provided within this version. Please see
e.g. Yu et al. (2009) for this means.

