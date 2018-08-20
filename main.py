'''
Anna helped me with my code 
Mitra helped me with my code 
'''

def main ():
  string1IO = open('Shot.txt','r')#choose the first protein sequence
  string1 = string1IO.read()
  string1=string1.replace(' ','').replace('\n','').strip()#strips blank spaces and new lines
  string1=string1[:300]#sets how long we want the sequence to be 
  string1IO.close()
 
  string2IO = open('hACF7.txt','r')#choose the second protein sequence 
  string2 = string2IO.read()
  string2=string2.replace(' ','').replace('\n','').strip()
  string2=string2[:300]
  string2IO.close()
  
  matrixName='BLOSUM62.txt'# choose what matrix you wanna use
  pam=getPAM(matrixName)#makes the scoring matrix(doesnt have to be pam)
  indel=-4#got this value from the book/rosalind, -8 is the indel value to use if using a PAM250 matrix, -4 if BLOSUM62
  
  print('string1 \n',string1) #just printing the inputs
  print('\n string2 \n',string2)
 
  finalScore,align1,align2,slimilarity,identity=localAlign(string1,string2,indel,pam)#gets all the statistics from the globalAlign function
  
 #summary stats go here. i got the round function off of stack overflow
  print('\nAlign Length:',max(len(align1),len(align2)),'\nScore:',finalScore,'\nSimilarity: ',round((slimilarity/max(len(align1),len(align2)))*100,2),'%','\nIdentity: ',round((identity/max(len(align1),len(align2)))*100,2),'%\n')
  makeAlign(align1,align2)#prints the alignment to the console
  return

def alignSymbol(align1,align2):#choose the symbols you want for your alignment here!
  alignMid=''#the middle alignment, where the symbols go 
  for i in range(len(max(align1,align2))):
    if align1[i]==align2[i]:
      alignMid=alignMid+align1[i] #'symbol' for exact match, but dont change this because if its exact match there doesnt need to be a symbol 
    elif align1[i]=='-':
      alignMid=alignMid+'_' #symbol for indel
    elif align2[i]=='-':
      alignMid=alignMid+'_' #symbol for indel
    else:
       alignMid=alignMid+'+' #symbol for similarity
  return alignMid
  
def makeAlign(align1,align2):#formats everything so nice and pretty into three row blocks like BLAST outputs
  alignMid=alignSymbol(align1,align2)
  maxLen=max(len(align1),len(align2))
  count=50#how long do you want the blocks?
  for i in range(0,maxLen,count):
    if i==0: #at the start start at 0 end at i+count. i learned of zfill on stack overflow.
      print(str(i).zfill(4),align1[i:(i+count)],str(i+count).zfill(4))
      print('    ',alignMid[i:(i+count)])
      print(str(i).zfill(4),align2[i:(i+count)],str(i+count).zfill(4))
      print()
    elif i+count < maxLen: #in the middle start at i+1 end at i+count
      print(str(i+1).zfill(4),align1[i:(i+count)],str(i+count).zfill(4))
      print('    ',alignMid[i:(i+count)])
      print(str(i+1).zfill(4),align2[i:(i+count)],str(i+count).zfill(4))
      print()
    elif i+count > maxLen: #at the end end at the end
      print(str(i+1).zfill(4),align1[i:(i+count)],str(maxLen).zfill(4))
      print('    ',alignMid[i:(i+count)])
      print(str(i+1).zfill(4),align2[i:(i+count)],str(maxLen).zfill(4))
      print()
  return

def getPAM(matrixName):

  '''
  Given to me by Anna
  Gets the PAM matrix by reading it from the PAM250.txt file.
  Input: nothing
  Output: dictionary of dictionaries: you can access scores like this:
  pam[’W’][’Y’] # score of a W/Y mismatch
  pam[’D’][’Y’] # score of a D/Y mismatch
  pam[’Y’][’Y’] # score of a Y/Y match
  '''
# get list of lines
  
  fin = open(matrixName)
  lines = fin.readlines()
  fin.close()
  pam = {}
  
  # go through each line from PAM250 file
  
  firstLine = True
  for line in lines:
  # row is now a list of strings
    row = line.strip().split()
  
    if firstLine:
    
    # if this is the first line, we have the headers.
    # store in the header list and set firstLine to false.
      header = row
     
      firstLine = False
  
  
    else:
    
    # this is a row of scores. The first element is the character.
      chari = row[0]
    # build a dictionary with (key,value) pairs, where the
    # keys are indexed from the header file. The values are ints.
      pam[chari] = {}
    
      for j in range(1,len(row)):
        charj = header[j-1]
        pam[chari][charj] = int(row[j])
  # printing dictionary; use header to order the elements

  for i in header:
    row = []
    for j in header:
      row = row + [pam[i][j]]

  
  print()
  return pam

def initializeTable(string1,string2): #taken from lab 8 
  table=[]
  for i in range(len(string1)+1):
      table=table+[[0]*(len(string2)+1)]
  return table
  
def m (c1,c2,pam):
    score= pam[c1][c2]
    #print('score: ',score)
    return score

def getSimilar (c1,c2,pam,slimilarity):
    
    score= pam[c1][c2] #if the PAM matrix score is positive, increase slimilarity counter by one
    if score >0:
      slimilarity +=1
    
    return slimilarity


def getIdentity (c1,c2,identity):
  #if the residues are identical, increase identity counter by 1
    if c1 == c2:
      identity +=1
    
    return identity


#this is all stolen from homework 6.2
def globalAlign(string1,string2,indel,pam):
    slimilarity=0
    identity=0
    align1 = ''
    align2 = ''
    table=initializeTable(string1,string2)
    backtrack=initializeTable(string1,string2)
    #print(table)
    for i in range(len(table)):
      for j in range(len(table[i])):
        holdingPen=[] #a place to take the max
        if i> 0:
          holdingPen=holdingPen+[table[i-1][j]+indel]#prev row
        if j >0:
          holdingPen=holdingPen+[table[i][j-1]+indel]#prev column
        if i > 0 and j > 0:
          holdingPen=holdingPen+[table[i-1][j-1]+m(string1[i-1],string2[j-1],pam)]#prev row and column
        if len(holdingPen) > 0:
          table[i][j]=max(holdingPen)
        #print('holdingPen:',holdingPen)
        #making the backtrack table now 
        if i == 0 and j == 0:
          backtrack[i][j]='*'
        if i >0 and table[i][j] == table[i-1][j]+indel:
          backtrack[i][j]='S'
        if j>0 and  table [i][j] == table[i][j-1]+indel:
          backtrack[i][j]='E'
        if i>0 and j>0 and  table [i][j] == table[i-1][j-1]+m(string1[i-1],string2[j-1],pam):
          backtrack[i][j]='D'
      #print(backtrack)
    #print(table)
    i = len(string1)
    j = len(string2)
    while i >0 or j > 0:
      if backtrack[i][j]=='S':#go up a row, add a dash to align2, the character to align1
        i += -1
        j += 0
        align2 = '-'+align2
        align1 = string1[i]+align1
      if backtrack[i][j]=='E':#go back a column, add a dash to align1, the character to align2
        i += 0
        j += -1
        align2 = string2[j]+align2
        align1 = '-'+align1
      if backtrack[i][j]=='D':#go back and up and add the characters to both aligns
        i += -1
        j += -1
        align2 = string2[j]+align2
        align1 = string1[i]+align1
        slimilarity=getSimilar(string1[i-1],string2[j-1],pam,slimilarity)
        identity=getIdentity(string1[i-1],string2[j-1],identity)
    finalScore = table[len(string1)][len(string2)]
   
    return finalScore,align1,align2,slimilarity,identity

'''
I attempted to get the local alignment working with my own code
first but i could not figure out where i went wrong so to test 
things I copied the whole of the model solution for local 
alignment to see if the error is in my attempt at integrating
the scoring matrix or in the my old code itself
'''

def localAlign(s1,s2,indel,pam):
  # initialize table with s1+1 rows and
  # s2+1 columns (remember letters span city 
  # blocks but we want nodes at intersections)
  identity=0
  similarity=0
  table = initializeTable(s1,s2)
  backtrack = initializeTable(s1,s2)
  
  for i in range(len(table)):
    for j in range(len(table[i])):
      if i==0 and j==0:
        table[i][j] = 0
        backtrack[i][j] = '*'
      if i>0 and j==0:
        table[i][j] = table[i-1][j]+indel
        backtrack[i][j] = 's'
      if i==0 and j>0:
        table[i][j] = table[i][j-1]+indel
        backtrack[i][j] = 'e'
      if i>0 and j>0: 
        table[i][j] = max([table[i-1][j-1]+m(s1[i-1],s2[j-1],pam)])
        if table[i][j] == table[i-1][j]+indel:
          backtrack[i][j] = 's'
        elif table[i][j] == table[i][j-1]+indel:
          backtrack[i][j] = 'e'
        else: # table[i][j] must equal table[i-1][j-1]+m]
          backtrack[i][j] = 'd'
      ## if table[i][j] < 0, then reset it.
      # this is the free taxi ride.
      if table[i][j]<0:
        table[i][j] = 0
        backtrack[i][j] = '*'

  # get max val  
  max_val = 0
  i = 0
  j = 0
  for this_i in range(len(table)):
    for this_j in range(len(table[this_i])):
      if table[this_i][this_j] > max_val:
       max_val = table[this_i][this_j]
       i = this_i
       j = this_j
  score = table[i][j]
  
  # now compute backtrack.
  align1 = ''
  align2 = ''
  while backtrack[i][j] != '*':
    print(backtrack[i][j])
    if backtrack[i][j] == 'd':
      align1 = s1[i-1]+align1
      align2 = s2[j-1]+align2
      i = i-1
      j = j-1
    elif backtrack[i][j] == 'e':
      align1 = '-'+align1
      align2 = s2[j-1]+align2
      j = j-1
    else: # backtrac[i][j] MUST be 's'
      align1 = s1[i-1]+align1
      align2 = '-'+align2
      i = i-1
    similarity=getSimilar(s1[i-1],s2[j-1],pam,similarity)
    identity=getIdentity(s1[i-1],s2[j-1],identity)



  return score,align1,align2,similarity,identity


'''
this was my attempt at working in local alignment using my own code
'''

def localAlignAttempt1(string1,string2,indel,pam):
  slimilarity=0
  identity=0
  align1 = ''
  align2 = ''
  table=initializeTable(string1,string2)
  backtrack=initializeTable(string1,string2)
  for i in range(len(table)):
    for j in range(len(table[i])):
      holdingPen=[0] #a place to take the max, plus taxi ride
      if i> 0:
        holdingPen=holdingPen+[table[i-1][j]+indel]#prev row
      if j >0:
        holdingPen=holdingPen+[table[i][j-1]+indel]#prev column
      if i > 0 and j > 0:
        holdingPen=holdingPen+[table[i-1][j-1]+m(string1[i-1],string2[j-1],pam)]#prev row and column
      if len(holdingPen) > 0:
        table[i][j]=max(holdingPen)
      #making the backtrack table now 
      if i == 0 and j == 0:
        backtrack[i][j]='*'
      if i >0 and table[i][j] == table[i-1][j]+indel:
        backtrack[i][j]='S'
      if j>0 and  table [i][j] == table[i][j-1]+indel:
        backtrack[i][j]='E'
      if i>0 and j>0 and  table [i][j] == table[i-1][j-1]+m(string1[i-1],string2[j-1],pam):
        backtrack[i][j]='D'
    print(backtrack)
  '''
  Hayden gave me this:
  '''
  i=0
  j=0
  maxVal=0
  for x in range(len(table)):
    for y in range(len(table[0])):
      #print(x,y,maxVal)
      if table[x][y]>maxVal:
        maxVal=table[x][y]
        i=x
        j=y 
  while backtrack[i][j] != '*':
    if backtrack[i][j]=='S':#go up a row, add a dash to align2, the character to align1
      i += -1
      j += 0
      align2 = '-'+align2
      align1 = string1[i]+align1
    if backtrack[i][j]=='E':#go back a column, add a dash to align1, the character to align2
      i += 0
      j += -1
      align2 = string2[j]+align2
      align1 = '-'+align1
    if backtrack[i][j]=='D':#go back and up and add the characters to both aligns
      i += -1
      j += -1
      align2 = string2[j]+align2
      align1 = string1[i]+align1
  finalScore = maxVal
  
  slimilarity=getSimilar(string1[i-1],string2[j-1],pam,slimilarity)
  identity=getIdentity(string1[i-1],string2[j-1],identity)
   
  
  return finalScore,align1,align2,similarity,identity


main()