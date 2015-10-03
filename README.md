# PromptSecondaryLTUnb

This is the repository for an attempt to create a prompt-secondary charm classifier which is
not correlated with the lifetime of the charm particle which you want to classify.

Directory structure : 

data   : contains the raw ntuples used for training

python : contains the code used to reduce ntuples into training format, eventually will contain
         scikit-learn based code for training classifiers

scripts : contains TMVA code used to train/test/create the classifier

classifiers : a repositry of trained classifiers with information about how they were trained
