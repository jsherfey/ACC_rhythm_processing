# dnsim-database-additions
- This is my personal store of homemade DNSIM mechanisms and models.
- As of now, this is Version 0.1, aka stuff still needs to have more descriptive comments added, and better explanation of everything.
- This is a more modern implementation of, basically, the thalamus component of
  - "Ching, S., Cimenser, A., Purdon, P. L., Brown, E. N., & Kopell, N. J. (2010). Thalamocortical model for a propofol-induced -rhythm associated with loss of consciousness. Proceedings of the National Academy of Sciences, 107(52), 22665-22670. http://doi.org/10.1073/pnas.1017069108"
- The two models currently in here, `TC-RE.mat` and `TC-RE-D96syns.mat` contain the basic TC-RE model with very few inputs (like from `itonic.txt`, `icosine.txt`, or `isquare.txt`) with merely different sets of synaptic strengths.
- This repo will receive more attention in the future.

- To use this:
  - Requires a working knowledge of DNSIM. Once you have that, you should be able to put the `*.mat` and `*.txt` files in your `dnsim/database` and go.
  - Or, to construct and run the model from the command line to do "industrial"-sized simulation like I do, run the example script `p1d1q023c1_D96_syn_ratios.m` via a simple `matlab -r p1d1q023c1_D96_syn_ratios` once your working directory is inside your `dnsim` install folder.

- If you find any typos, **PLEASE** let me know!!!
