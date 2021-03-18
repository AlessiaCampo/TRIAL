
'''
This program allows to compute the best local alignment between two aminoacid sequences respect to the scoring system that includes a substitution matrix and a gap penalty,
based on the Smith-Waterman alogrithm.

1.INITIALIZE AND FILL THE MATRICES F AND P:
Define a function fill_and_score() that takes in input:
-the two aminoacid sequences (seq1 and seq2);
-the gap penalty (p)
-the substitution matrix (SM).
        Initialize the scoring matrix F and the derivation matrix P  //list comprhension
        One-to one comparison of residues in seq1 and seq2
        Score pairs of residues using the SM and gap penalty (p) (adding gaps or substitution scores)
        Associate scores with left, up and diagonal derivations   //zip built-in function
        Fill F with the best scores among the three possible outcomes (gap from the left, gap from up, substitution score from the diagonal) for each pair of residues  //max buil-in function
        Fill P with the derivations left, up and diagonal respect to the score assigned in the F matrix
        Add the condition that negative scores will be equal to zero in the scoring matrix
    Return F and P filled

2. TRACEBACK AND FINAL OUTPUT
Define a function best_local() that takes in input:
-the matrices F and P
-the two sequences seq1 and seq2
        Select the highest score (max) and its position in the matrix F
        Start the traceback from the max according to the path of the P matrix until zero is found // use indexes to move trhough the matrix
    Return strings containing the residues aligned and the final score (max)
'''
#IMPORT THE SUBSTITUTION MATRIX AND THE AMINOACID SEQUENCES
import input_data

SM=input_data.BLOSUM52 #import the substitution matrix
#print(SM)
seq1=input_data.seq1  #import the first aminoacid sequence
#print(seq1, len(seq1))
seq2=input_data.seq2  #import the second aminoacid sequence
#print(seq2,len(seq2))

pen=-2
def fill_and_score(pen, seq1, seq2, SM):
#INITIALIZE THE MATRICES F AND P
    F=[[0]*(len(seq1)+1) for j in range(len(seq2)+1)]   #initialize the F matrix as a list of lists where the rows are defined by the length of seq2 +1 and the number of columns by the length of seq1+1
    #print(F)
    P=[["0"]*(len(seq1)+1) for i in range(len(seq2)+1)] #initialize the derivation matrix P as list of lists as done for the F matrix
    #print(P)
#SCORING:FILL THE MATRICES F AND P
    for col in range(1,len(seq1)+1):  #iterate over the length of seq1 starting from the first position (the first column is not filled)
        #print("col",col)
        for row in range(1,len(seq2)+1):  #iterate over the length of seq2 starting from the first position (the first row is not filled)
            #print("row", row)
            pair=seq1[col-1]+seq2[row-1] #create pairs of residues by concatenating the letters from seq1 and seq2
            diag=SM[pair]+F[row-1][col-1] #define pairs as keys of the dictionary SM. The values associated to each keys of the dictionary are the match scores added to the score derived from the diagonal position of the F matrix
            left=F[row][col-1]+pen #the score assigned to pairs  computed by the addition of the gap penalty (mismatch) to the score of the left position in the F matrix
            up=F[row-1][col]+pen #the scores assigned to pairs computed by the addition of the gap penalty (mismatch) to the score of the up position in the F matrix
            pair_scores=list(zip((diag,left,up),("d","l","u"))) #for each pair of residues the zip function allow to associate each type of score defined above to the strings of the derivation respectively; the function list converts these groups of associated values into a lists
            #print(pair_scores)
            F[row][col],P[row][col]=max(pair_scores) #using the max built-in function, the F matrix is filled with the maximum scores selected among the three different ones and the P matrix is filled with the corresponding string of derivation
            #print(max(ass))
            if F[row][col]<0: #this condition is used to fill the matrix F with zero when the score selected is negative since negatives score are not accepted in the local alignment
                F[row][col]=0
                P[row][col]="0"

    return(F,P)


F,P=fill_and_score(pen,seq1,seq2,SM)
#print(F)
#print(P)

def best_local(F,P, seq1, seq2):
#TRACEBACK1: DEFINE THE HIGHEST SCORE POSITION
    max=F[0][0]                     #define the first position of the F matrix as the maximum value to set a starting point and find the highest score considering all the positions
    for col in range(1,len(seq1)+1): #iterate over the length of the seq1+1
        for row in range(1,len(seq2)+1): #iterate over the length of the seq2+1
            if F[row][col]>max: #the condition is used to set the highest score in the iteration. At each iteration the value of max is updated untill the highest value is found
                c=col #assign to a new variable c the coordinate col of the maximum value
                r=row #assign to a new variable r the coordinate row of the maximum value //this is done to save the position of the highest score
                max=F[r][c] #when the highest score is found max take that value in the corresponding position over the matrix F. It is the score of the best local alignment

#TRACEBACK2: FOLLOW THE PATH BASED ON DERIVATIONS
    s1="" #define two empty strings s1 and s2 in which insert the residues aligned
    s2=""
    while F[r][c]!=0: # iterate over the matrices positions until a zero in the F matrix is found, starting from the coordinates saved previously as c and r of the highest value
        if P[r][c]=="l": #1st condition: when "l" is found in the P matrix:
            s1=s1+seq1[c-1] #the string s1 is filled with the corresponding letter in the seq1 in the current position(matrix F has an extra position respect to the sequences indexes)
            s2=s2+"-" #simultaneously in the string s2 a gap is inserted
            c=c-1 #specify the position for the next iterations: one columns before following the left derivation is the position to be iterated in the next cycle  //traceback
        elif P[r][c]=="u": #2nd condition: when "u" is found in the P matrix:
            s1=s1+"-" #a gap is inserted in the first sequence s1
            s2=s2+seq2[r-1] #simultanously the  string s2 is filled with the corresponding letter in the seq2 in the current position (matrix has an extra position )
            r=r-1  #specifiy the position for the next iteration: one row before following the up derivation is the position to beiterated in the next cycle //traceback
        elif P[r][c]=="d": #3th condition: when "d" is found in the P matrix:
            s1=s1+seq1[c-1] #the first string is filled with the letter in the corresponding position in seq1 (matrix has an extra position)
            s2=s2+seq2[r-1] #the second string s2 is filled with the matching letter in the position in seq2 (matrix has an extra position)
            c=c-1  #one column and row before following the diagonal derivation is the position to be iterated in the next cycle //traceback
            r=r-1


    return(s1[::-1],s2[::-1], max) #print the string in reverse to visualize them in the correct order


#(best_local(F,P, seq1,seq2))
s1,s2,max=best_local(F,P,seq1, seq2)
#FINAL OUTPUT
print ("BEST LOCAL ALIGNMENT: ")
print( "score:" , max)
print(s1)
print(s2)
