function ensemble=savePhy()

[cids,cgs]=readClusterGroupsCSV('cluster_group.tsv');
goodUnits=cids((cgs==2));%good units
clusters=readNPY('spike_clusters.npy');
ts=readNPY('spike_times.npy');

ensemble=cell(length(goodUnits),3);
for i=1:length(goodUnits)
	ensemble{i,3}=ts(clusters==goodUnits(i))';
end

return
