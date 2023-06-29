# Attune_Cruise_Process
Scripts for processing attune flow cytometry data from Sosik Lab cruises, underway and discrete sampling


Below is the protocol for processing discrete samples. 

# Attune Discrete Sample Protocol 

The main script for processing a cruise of discrete samples is "Process_Preserved_AllSteps.m"

To run this script, you will need to adjust the top for a specific cruise: 
    basepath = where to look for preserved fcm sample fcs fles
    cruisename 
    restpath = where to look for CTD data as metadata for niskin files 
    elogpath = where to look for metadata associated with underway discrete files if there are any

    you can also adjust the flags Step1 and Step5only if you don't need to redo file set up (step 1) or if you only want to redo size calibration (step 5)
    if you are processing a cruise from the beginning, set Step1 ==1 and Step5only == 0. 

    hierarchical_gates can be set to either 'True' or 'False'. This has to do with whether gates were drawn within other gates. 
        early cruises did not have hierarchical gating, but OTZ cruises that Alexi gated do, and subsequent cruises usually will have hierarchical gates. 
        As of 6/29/23, code will only look for a single parent gate. No double parents. Grandparents are fine. 


Then run the script. Fingers crossed it all goes smoothly. 



