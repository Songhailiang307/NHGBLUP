# -*- encoding: utf-8 -*-
import sys

up_triangle = open(sys.argv[1]+'_dmu', 'w')
down_triangle = open(sys.argv[1]+'_asreml', 'w')

ID={}
e=1
with open(sys.argv[2]) as id:
	for i in id:
		f=i.split()
		ID[str(e)]=f[0]
		e+=1

with open(sys.argv[1]) as ff:
	for n,line in enumerate(ff):
		line = line.split()
		for i in range(len(line)):
			if i <= n:
				if sys.argv[3]=='Gma':
					if i == n:
						val = float(line[i]) + 0.001
					else:
						val = line[i]
					down_triangle.write('%s\t%s\t%s\n' % (ID[str(n+1)], ID[str(i+1)], str(float(val))))
				elif sys.argv[3]=='Ginv':
					val = line[i]
					down_triangle.write('%s\t%s\t%s\n' % (ID[str(n+1)], ID[str(i+1)], str(float(val))))
				else:
					print ("error!!")
					break
			if i >= n:
				if sys.argv[3]=='Gma':
					if i == n:
						val = float(line[i])+0.001
					else:
						val = line[i]
					up_triangle.write('%s\t%s\t%s\n' % (ID[str(n+1)], ID[str(i+1)], str(float(val))))
				elif sys.argv[3]=='Ginv':
					val = line[i]
					up_triangle.write('%s\t%s\t%s\n' % (ID[str(n+1)], ID[str(i+1)], str(float(val))))
				else:
					print ("error!!"	)
					break
					
up_triangle.close()
down_triangle.close()
                