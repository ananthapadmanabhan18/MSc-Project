# Running the Analysis

<ol>
    <li> <b> Run the MCInfo:</b> Edit the Paths (simFolder,outputMCIFolder etc) in the MCInfo.C and then run the runMCInfo.sh script by <code> ./runMCInfo n1 n2 ini_event_no final_event_no</code> </li>
    <li> <b> Run the Clusteriser:</b> Edit the Paths (simFolder,outputMCIFolder etc) in the Analysis.C and then run the runAnalysis.sh script by <code> ./runclustering n1 n2 ini_event_no final_event_no</code></li>
    <li> <b> Run the Analysis:</b> Edit the Paths (simFolder,outputMCIFolder and also the Location of <code>MCInfoXXXX.root</code> and <code>clusterXXX.root</code>) in the Analysis.C and then run the runAnalysis.sh script by <code> ./runAnalysis n1 n2 ini_event_no final_event_no</code></li>
</ol>