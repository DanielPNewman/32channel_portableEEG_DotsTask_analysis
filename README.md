# 32channel portable EEG Dot Task analysis
Load process and analyse human 32 channel EEG data collect using BrainProducts’ portable LiveAmp EEG system. This data was collected while participants performed a variant of the random dot motion task (e.g. Britten, Shadlen, Newsome, & Movshon, 1992; Kelly & O’Connell, 2013; Loughnane et al., 2016; Newsome, Britten, & Movshon, 1989)

Borrows a lot of code from this repo https://github.com/gerontium/big_dots

### Run the scripts in this order:

1. getChannelVars.m
2. Check4badchans.m
3. runafew.m
4. ERP_analysis.m