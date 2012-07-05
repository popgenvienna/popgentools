import random
import unittest
import CMH
import RCMH




class TestSyncReader(unittest.TestCase):
	"""
	Thoroughly test the reade for sync files
	Create an artificial sync file as a string
	open a filehandle to it and check if the output is fine
	""" 
	def setUp(self):
		import StringIO
		self.file=StringIO.StringIO(
"""
2L	1	C	1:2:3:4:5:6	0:6:8:0:0:1
2L	2	C	0:20:4:0:0:1	1:0:4:0:0:4
2L	3	C	4:0:0:0:0:0	4:0:0:0:0:0
""")
		
	def test_reader(self):
		file=self.file
		sr=CMH.SyncReader(file)
		l=sr.next()
		self.assertEqual(l.chr,"2L")
		self.assertEqual(l.popcount,2)
		self.assertEqual(l.pos,1)
		self.assertEqual(l.refc,"C")
		p1=l.populations[0]
		self.assertEqual(p1.A,1)
		self.assertEqual(p1.T,2)
		self.assertEqual(p1.C,3)
		self.assertEqual(p1.G,4)
		self.assertEqual(p1.N,5)
		self.assertEqual(p1.deletion,6)
		p2=l.populations[1]
		self.assertEqual(p2.A,0)
		self.assertEqual(p2.T,6)
		l=sr.next()
		self.assertEqual(l.chr,"2L")
		self.assertEqual(l.pos,2)
		self.assertEqual(l.refc,"C")
		l=sr.next()
		self.assertEqual(l.pos,3)
		self.assertRaises(StopIteration,sr.next)


class TestSinglePop(unittest.TestCase):
	def setUp(self):
		self.pop = CMH.SinglePop(1,2,3,4,5,6)
		
	def test_properties(self):
		p=self.pop
		self.assertEqual(p.A,1)
		self.assertEqual(p.T,2)
		self.assertEqual(p.C,3)
		self.assertEqual(p.G,4)
		self.assertEqual(p.N,5)
		self.assertEqual(p.deletion,6)
		self.assertEqual(p.cov,10)
		self.assertEqual(p.totcov,21)
		
		
		

class TestRCMHUtil(unittest.TestCase):
	def setUp(self):
		self.p1=[CMH.SinglePop(3,1,4,0,0,0),CMH.SinglePop(3,4,0,4,8,12)]
		self.p2=[CMH.SinglePop(3,4,0,0,0,0),CMH.SinglePop(3,4,0,0,0,0),CMH.SinglePop(0,0,6,7,0,0),CMH.SinglePop(2,1,1,0,0,0)]
		self.majorcount1=RCMH.Utility.getMajorAlleleCount(self.p1)[0]
		self.majorcount2=RCMH.Utility.getMajorAlleleCount(self.p2)[0]
		
	
	def test_getMajorAlleleCount(self):
		mc=self.majorcount1
		self.assertEqual(mc[0][0],3)
		self.assertEqual(mc[1][0],3)
		self.assertEqual(mc[0][1],1)
		self.assertEqual(mc[1][1],4)
		mc=self.majorcount2
		self.assertEqual(mc[0][1],3)
		self.assertEqual(mc[1][1],3)
		self.assertEqual(mc[2][1],0)
		self.assertEqual(mc[3][1],2)
		self.assertEqual(mc[0][0],4)
		self.assertEqual(mc[1][0],4)
		self.assertEqual(mc[2][0],0)
		self.assertEqual(mc[3][0],1)

	def test_issnp(self):
		self.assertTrue(RCMH.Utility.issnp(self.majorcount1,5))
		self.assertFalse(RCMH.Utility.issnp(self.majorcount1,6))
		self.assertTrue(RCMH.Utility.issnp(self.majorcount2,8))
		self.assertFalse(RCMH.Utility.issnp(self.majorcount2,9))
	
	def test_mantelhaen(self):
		
		# no change
		pvalue=RCMH.Utility.getCMHPvalue(((2,2),(4,4),(3,3),(5,5)))
		self.assertTrue(abs(pvalue-1)<0.0001)
		
		# trends with variing numbers of replicates
		pvalue=RCMH.Utility.getCMHPvalue(((10,0),(10,5),(40,0),(40,20)))
		self.assertTrue(abs(pvalue-1.7e-5)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((10,0),(10,5),(10,0),(10,5)))
		self.assertTrue(abs(pvalue-0.01333)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((10,0),(10,5),(10,0),(10,5),(10,0),(10,5)))
		self.assertTrue(abs(pvalue-0.001496)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((10,0),(10,5),(10,0),(10,5),(10,0),(10,5),(10,0),(10,5)))
		self.assertTrue(abs(pvalue-0.0001768)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((20,0),(20,10),(20,0),(20,10)))
		self.assertTrue(abs(pvalue-0.0001513)<0.0001)
		
		# lets try opposing trends
		pvalue=RCMH.Utility.getCMHPvalue(((20,20),(20,40),(20,20),(40,20)))
		self.assertTrue(abs(pvalue-1)<0.0001)
		
		# shrinking values
		pvalue=RCMH.Utility.getCMHPvalue(((20,20),(20,10),(20,20),(20,10)))
		self.assertTrue(abs(pvalue-0.07401)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((20,20),(20,10),(20,20),(20,10),(20,20),(20,10)))
		self.assertTrue(abs(pvalue-0.02394)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((20,20),(20,0),(20,20),(20,0)))
		self.assertTrue(abs(pvalue-1.716e-07)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((20,20),(20,0),(20,20),(20,0),(20,20),(20,0)))
		self.assertTrue(abs(pvalue-8.828e-11)<0.0001)
		pvalue=RCMH.Utility.getCMHPvalue(((10,90),(100,0),(10,90),(100,0)))
		self.assertTrue(abs(pvalue-5.2373623318431986e-72)<0.0001)

class TestRCMHSimulate(unittest.TestCase):

	
	def test_simulateDriftSingleSample(self):
		s=RCMH.SimulateDrift.singleSampleDrift(1000,(50,50),(10,20,30,40),(30,40,60,70))
		
		s=RCMH.SimulateDrift.singleSampleDrift(1000,(100,0),(10,20,30,40),(30,30,30,30))
		
		s=RCMH.SimulateDrift.singleSampleDrift(1000,(0,100),(10,20,30,40),(30,30,30,30))
		
	def test_createMatrix(self):
		m=RCMH.SimulateDrift._createMatrix((1,2,3,4,5,6),3)
		self.assertEqual(m[0][0],1)
		self.assertEqual(m[0][1],4)
	
	def test_matrixToList(self):
		m=RCMH.SimulateDrift._matrixToList(((10,40),(20,50),(30,60)))
		self.assertEqual(m[0],10)
		self.assertEqual(m[5],60)


		
		
	def test_multisampledrift(self):
		s=RCMH.SimulateDrift.multiSampleDrift(1000,((20,20),(0,20),(20,0)),(20,20,20,40,40,40),(10,20,30,40,50,60))

		#1000,((20,20),(0,20),(20,0)),(20,20,20,40,40,40),(10,20,30,40,50,60) -> [(8, 2), (0, 20), (30, 0), (30, 10), (0, 50), (60, 0)]
		bla=1



if __name__ == '__main__':
	unittest.main()