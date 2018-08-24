def motif():
    x=input('Please input your sequence motif capitalized that you are looking for (X indicates any nucleotide): ')
    y=input('Please input your sequence (capitalized and no spaces) that you would like to search: ')
    indlst=[]#list of indexes of the nucleotides present in motif
    letlst=[]#list of nucleotides for the motif
    count=0
    for let in x:#loop to get the index and nucleotides for the motif
        if let=='A' or let=='G' or let=='C' or let=='T' or let=='U':
            indlst.append(count)
            letlst.append(let)
            count=count+1
        else:
            count=count+1
    seqlst=[]
    for i in range(((len(y)-len(x)))):#moving window for motif search--this window is same size as motif
        c=len(x)+i
        seq=y[i:c]
        seqlst.append(seq)#add every sequence of the moving window to a list
    finalseqlst=[]
    for letter in seqlst: #get the relevant nucleotides from the sequence, the ones with the same index, for every item
        newstr=''
        for num in indlst:
            newstr=newstr+letter[num]
        finalseqlst.append(newstr)
    countlst=[]
    newcount=0
    newnewcount=0
    z=''.join(letlst)
    for a in finalseqlst:#compare the motif and the and sequence
        if a==z:
            countlst.append(newnewcount)
            newcount=newcount+1
            newnewcount=newnewcount+1
        else:
            newnewcount=newnewcount+1
    print(newcount)
    returnseqlst=[]
    for b in countlst:#get the original sequence, not just the indexed sequence
        returnseqlst.append(seqlst[b])
    print(returnseqlst)

    #to translate the sequence to the negative strand and search
    negstrndstr=''
    for nuc in y:
        if nuc=='A':
            negstrndstr=negstrndstr+'T'
        if nuc=='G':
            negstrndstr=negstrndstr+'C'
        if nuc=='C':
            negstrndstr=negstrndstr+'G'
        if nuc=='T':
            negstrndstr=negstrndstr+'A'
        if nuc=='U':
            negstrndstr=negstrndstr+'A'
    print(negstrndstr)
    negseqlst=[]
    for i in range(((len(negstrndstr)-len(x)))):
        f=len(x)+i
        negseq=negstrndstr[i:f]
        negseqlst.append(negseq)
    finalnegseqlst=[]
    for nucls in negseqlst:
        negnewstr=''
        for negnum in indlst:
            negnewstr=negnewstr+nucls[negnum]
        finalnegseqlst.append(negnewstr)
    negcountlst=[]
    negnewcount=0
    negnewnewcount=0
    for h in finalnegseqlst:
        if h==z:
            negcountlst.append(negnewnewcount)
            negnewcount=negnewcount+1
            negnewnewcount=negnewnewcount+1
        else:
            negnewnewcount=negnewnewcount+1
    print(negnewcount)
    negreturnseqlst=[]
    for e in negcountlst:
        negreturnseqlst.append(negseqlst[e])
    print(negreturnseqlst)

        
    
motif()
        
    
    
                   
                       

                       
                   
            
        
    
                   
                   
    
