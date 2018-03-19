# clas-phi
This repo holds tools and documentation for the study of phi (1020) from E1-F.  

### Setup 
- First you should have ROOT installed 
- Clone this repo 
- cd data/ && wget https://userweb.jlab.org/~dmriser/lfs/events.root

### Structure of Data
This data has been skimmed with (W > 2, Q2 > 1).  It also has loose skims on the mass of the missing state in the reactions to reduce the size of the file.  These can be adjusted at request.

- helicity - beam helicity 
- topology - 3 is fully exclusive, 0 means no negative kaon was detected and the TLorentzVector for kaon_neg is the missing state X 
- alpha - the confidence level for each particle id, in the case of no particle it is set to 0.  This has been skimmed to be above 1e-3 for the detected particles
- dist - These distances are used to check systematics, and represent the distance away from nominal cut values which are -1 and 1 in the rescaled plot.  
- electron, proton, kaon_pos, kaon_neg - These are TLorentzVectors that have been identified by my other analysis routines.
