#!/usr/bin/python

import sys

names = ['1h7m','1ihb','1io2','1iro','1jiw','1kfw','3run','1lhm','1mjc','3snf','2qmt',
'1rbp','1rro','1rwy','3ua7','1u06','1sup','1pml','1urp','4n0k','1yea','1yu5','1zym','2hbb','2hip','2rn2',
'2trx','3d2c','3sil','451c','4blm', '1a23','3a9j','1ag2','1aj3','3wc8','1am7','1arr','4hil','1a2p','1bpi','1a6m','1c2r','1cdc',
'1clw','3pf4','1cyo','1el1','1rcf','3b0k','1hue','1ig5','2nvh','1k9q','2nwd','1mbg','1mgr','1otr','2qmt',
'1poh','3ne4','1lni','1rn1','2ijk','1kf5','1ra9','4wor','1gvp','1wq5','1hb6','1rg8','1aky','2ci2','2hpr',
'3fa0','2cle','2hds','1hwg','3q27','3ssi','2vb1','1gvp','1wq5','1hb6','1rg8','1aky','2ci2','2hpr',
'3fa0','2cle','2hds','1hwg','3q27','3ssi','2vb1','5cro','4jzk','2bt6','1azp','1c52','1c9o','1chk','1cyo',
'1v0l','1f6r','4wf5','1el1','3b0k','1ig5','1k9q','1mbg','1otr','2qmt','3ne4','2ijk','1kf5','2ci2','1hwg','2vb1',
'1azp','3run','3snf','2qmt','1rro','1rwy','1pml','2vgl']



for name in names:
	f = open(name+'NEW.txt','r')
	out = open(name+'_NEW.txt','w')
	for line in f.readlines():
		out.write(line[0])
		for c in line:
			if(c.isdigit()):
				out.write(c)
		out.write('\n')
	f.close()
	out.close()
