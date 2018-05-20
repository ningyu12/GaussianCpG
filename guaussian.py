import sys

delta =  4.9;
e = 2.71828182846
A = 256;
#k = float(1/1.3)*float(1/float(delta*2.506628253));
k = float(1/float(delta*2.506628253));
#k = float(1/1.29);
#k = 1;


print 'delta ---->',delta, '\tA=', A, '\tk=',k;
for x in range(40):
	print x, '\t', float(e**(float(-x*x)/float(2*delta*delta)) ), '\t',(float(k )* float(A)*float(e**(float(-x*x)/float(2*delta*delta)) )), '\t',round(float(k )* float(A)*float(e**(float(-x*x)/float(2*delta*delta)) ));


gsfilename = "gaussian.txt";
gsfile = open(gsfilename, 'w+');
for x in range(40):
	if round(float(k )* float(A)*float(e**(float(-x*x)/float(2*delta*delta)) )) > 0:
		gsfile.write("{}\n".format(int(round(float(k )* float(A)*float(e**(float(-x*x)/float(2*delta*delta)))))));
	else:
		break;

gsfile.close();
