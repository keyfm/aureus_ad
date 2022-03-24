def read_samplesCSV(spls):
	# reads in samples.csv file, format: Batch,Lane,Barcode,Sample,Alignments,ProviderName,Patient
	hdr_check = ['Batch', 'Lane', 'Barcode', 'Sample', 'Alignments', 'ProviderName', 'Patient']
	switch = "on"
	file = open(spls, 'r')
	list_path = []
	list_splID = []
	list_providerNames = []
	list_refG = []
	list_patient = []
	for line in file:
		line = line.strip('\n').split(',')
		# Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
		if switch == "on":
			if (line == hdr_check):
				print("Passed CSV header check")
			else:
				Warning("CSV did NOT pass header check! Code continues, but first line ignored")
			switch = "off"
			continue
		# build lists
		list_path.append(line[0])
		list_splID.append(line[3])
		list_refG.append(line[4])
		list_providerNames.append(line[5])
		list_patient.append(line[6])
	return [list_path,list_splID,list_refG,list_providerNames,set(list_patient)] # set(list_patient) provides only unique subject IDs


def read_samples_caseCSV(spls):
	# reads in samples_case.csv file, format: Path,Sample,ReferenceGenome,Outgroup
	hdr_check = ['Path','Sample','ReferenceGenome','Outgroup']
	switch = "on"
	file = open(spls, 'r')
	list_path = []
	list_splID = []
	list_refG = []
	list_outgroup = []
	for line in file:
		line = line.strip('\n').split(',')
		# Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
		if switch == "on":
			if (line == hdr_check):
				print("Passed CSV header check")
			else:
				Warning("CSV did NOT pass header check! Code continues, but first line ignored")
			switch = "off"
			continue
		# build lists
		list_path.append(line[0])
		list_splID.append(line[1])
		list_refG.append(line[2])
		list_outgroup.append(line[3])
	return [list_path,list_splID,list_refG,list_outgroup]

