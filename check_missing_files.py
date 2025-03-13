import os
import subprocess as sp

def main():
	directory = './'
	diroutput = './output/'
		
	output_index = []
	for filename in os.listdir(diroutput):
		if filename.startswith('output') and filename.endswith('root'):
			index_str = filename[-9:-5]			
			output_index.append(int(index_str))
	
	missing_index = []
	for filename in os.listdir(directory):
		if filename.startswith('sched') and filename.endswith('csh'):
			index_str = filename[filename.find('_')+1:-4]
			if int(index_str) not in output_index:
				# print(int(index_str))
				missing_index.append(int(index_str))

	print(len(missing_index))	
	missing_index_str = ''
	for index in missing_index:
		missing_index_str += f'{index},'
	missing_index_str = missing_index_str[:-1]

	# write condor file instead
	# get jobid, which should be in the format of 'sched5FB88B1A8CDB427602BCD52B2C699D47'
	jobid = filename.split('_')[0]	
	# gather pwd
	pwd = os.environ['PWD']
	with open('resubmit.condor', 'w') as f:
		for index in missing_index:
			f.write(f"""Universe         = vanilla
Notification     = never
Executable       = /bin/env
Arguments        = "singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif {pwd}/sched{jobid}_{index}.csh"
Output           = {pwd}/log/script_{index}.out
Log              = /tmp/maxwoo/{jobid}.condor.log
Initialdir       = {pwd}
kill_sig         = SIGINT
# PeriodicRemove   = (NumJobStarts >=1 && JobStatus==1) || (JobStatus == 2 && (CurrentTime - JobCurrentStartDate > (54000)) && ((RemoteUserCpu+RemoteSysCpu)/(CurrentTime-JobCurrentStartDate) < 0.10)) || (((CurrentTime - EnteredCurrentStatus) > (2*24*3600)) && JobStatus == 5) || (JobRunCount >= 1 && JobStatus == 1)
Priority         = +10
Queue


""")
			
	# with open('resubmit.sh', 'w') as f:
	# 	f.write('#!/bin/bash\n')
	# 	f.write(f'star-submit-beta -r {missing_index_str} *.session.xml')	


if __name__ == '__main__':
	main()
