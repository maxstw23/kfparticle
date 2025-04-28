import os
import subprocess as sp

def main():
	directory = './'
	diroutput = './output/'
		
	output_index = []
	jobid = ''
	for filename in os.listdir(diroutput):
		if filename.startswith('output_merged') and filename.endswith('root'):
			index_str = filename[filename.rfind('_')+1:-5]	
			output_index.append(int(index_str))
			jobid = filename[filename[:filename.rfind('_')].rfind('_')+1:filename.rfind('_')]
			
	missing_index = []
	# if unmerged_files2 exist in diroutput, check scripts there, otherwise check in unmerged_files1
	if os.path.exists(diroutput + 'unmerged_files2'):
		directory = diroutput + 'unmerged_files2/haddinfo/location/'
	elif os.path.exists(diroutput + 'unmerged_files1'):
		directory = diroutput + 'unmerged_files1/haddinfo/location/'
	else:
		print("No unmerged_files directory found.")
		return
    
	for filename in os.listdir(directory):
		if filename.startswith('hadd') and filename.endswith('csh'):
			index_str = filename[filename.find('_')+1:-4]
			if int(index_str) not in output_index:
				# print(int(index_str))
				missing_index.append(int(index_str))
	
	print(len(missing_index))	
	missing_index_str = ''
	for index in missing_index:
		missing_index_str += f'{index},'
	missing_index_str = missing_index_str[:-1]
	with open('resubmit_myhadd.sh', 'w') as f:
		f.write('#!/bin/bash\n')
		f.write(f'star-submit -r {missing_index_str} hadd{jobid}.session.xml')	


if __name__ == '__main__':
	main()
