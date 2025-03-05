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
	with open('resubmit.sh', 'w') as f:
		f.write('#!/bin/bash\n')
		f.write(f'star-submit-beta -r {missing_index_str} *.session.xml')	


if __name__ == '__main__':
	main()
