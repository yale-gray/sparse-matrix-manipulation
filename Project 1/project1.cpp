#include <iostream>
using namespace std;

class SparseRow{
    protected:
        int row;
        int col;
        int val;
    public:
        //default constructor
        SparseRow();
        //parametrized instructor
        SparseRow(int row, int col, int val);
        //display function
        void display();
        //overloads the output operator
        friend ostream& operator<<(ostream& s, const SparseRow);
        //getter functions
        int getRow() const;
        int getCol() const;
        int getVal() const;
        //setter functions
        void setRow(int r);
        void setCol(int c);
        void setVal(int v);
}; 

class IncompatibleMultiplicationException : public exception{
    public:
        const char* what() const throw(){
            //state when the two matrices are unable to by multiplied
            return "Matrix multiplication is not possible";
        }
};

class IncompatibleAdditionException : public exception{
    public:
        const char* what() const throw(){
            //state when the two matrices are unable to be added
            return "Matrix addition is not possible";
        }
};

SparseRow::SparseRow(){
    row = -1;
    col = -1;
    val = 0;
}
SparseRow::SparseRow(int row, int col, int val){
    this->row = row;
    this->col = col;
    this->val = val;
}
void SparseRow::display(){
    cout << row << ", " << col << ", " << val << endl;
}
ostream& operator<<(ostream& s, const SparseRow sr){
    s << sr.row << ", " << sr.col << ", " << sr.val << endl;
    return s;
}
//getters
int SparseRow::getRow() const{
    return row;
}
int SparseRow::getCol() const{
    return col;
}
int SparseRow::getVal() const{
    return val;
}
//setters
void SparseRow::setRow(int r){
    row = r;
}
void SparseRow::setCol(int c){
    col = c;
}
void SparseRow::setVal(int v){
    val = v;
}

class SparseMatrix {
    protected:
        int noRows; //Number of rows of the original matrix
        int noCols; //Number of columns of the original matrix
        int commonValue; //read from input
        int noNonSparseValues;
        SparseRow* myMatrix; //Array
        int matrixIndex;
        int** nonSparseMatrix; //Array of arrays in matrix form
    public:
        SparseMatrix ();
        SparseMatrix (int n, int m, int cv, int noNSV);
        SparseMatrix* Transpose (); //Matrix Transpose
        SparseMatrix* Multiply (SparseMatrix& M); //Used to multiple two matrices together
        SparseMatrix* Add (SparseMatrix& M); //Used to add two matrices together
        friend ostream& operator<< (ostream& s, const SparseMatrix& sm);
        void displayMatrix (); //Display the matrix in its original format
        void toMatrix(); //Helper method to convert to normal matrix form
        void readFile(SparseMatrix& M); //Helper method to read in values from a file
        //Getter functions
        int getNoRows() const;
        int getNoCols() const;
        int getCommonValue() const;
        int getNoNonSparseValues() const;
        SparseRow* getMyMatrix() const;
        int getVal(int r, int c) const;
        //Setter functions
        void setNoRows(int r);
        void setNoCols(int c);
        void setCommonValue(int cv);
        void setNoNonSparseValues(int n);
};

SparseMatrix::SparseMatrix(){
    noRows = -1;
    noCols = -1;
    commonValue = -1;
    noNonSparseValues = -1;
    myMatrix = NULL;
}
SparseMatrix::SparseMatrix(int n, int m, int cv, int noNSV){
    noRows = n;
    noCols = m;
    commonValue = cv;
    noNonSparseValues = noNSV;
    myMatrix = new SparseRow[noNonSparseValues];
}
SparseMatrix* SparseMatrix::Transpose(){
    //Create a new SparseMatrix object to store the transposed matrix
    SparseMatrix* sm = new SparseMatrix(noCols, noRows, commonValue, noNonSparseValues);
    for(int i = 0; i < noNonSparseValues; ++i){
        //Create a new SparseRow object to store the transposed row in the new matrix
        SparseRow* temp = new SparseRow(myMatrix[i].getCol(), myMatrix[i].getRow(), myMatrix[i].getVal());
        sm->myMatrix[i] = *temp;
    }
    return sm;
}
SparseMatrix* SparseMatrix::Multiply(SparseMatrix& M){
    //check that the matrices are compatible for multiplication
    if(noCols != M.getNoRows()){
        throw IncompatibleMultiplicationException();
    }
    //counter for the number of non-sparse values in the result matrix
    int resultNonSparseValues = 0;
    for(int i = 0; i < noNonSparseValues; ++i){
        for(int j = 0; j < M.getNoNonSparseValues(); ++j){
            if(myMatrix[i].getCol() == M.myMatrix[j].getRow()){
                ++resultNonSparseValues;
            }
        }

    }
    
    //create the result matrix
    SparseMatrix* result = new SparseMatrix(noRows, M.getNoCols(), commonValue, resultNonSparseValues);
    //fill the result matrix
    int resultIndex = 0;
    for(int i = 0; i < noNonSparseValues; ++i){
        for(int j = 0; j < M.getNoNonSparseValues(); ++j){
            if(myMatrix[i].getCol() == M.myMatrix[j].getRow()){
                result->myMatrix[resultIndex].setRow(myMatrix[i].getRow());
                result->myMatrix[resultIndex].setCol(M.myMatrix[j].getCol());
                result->myMatrix[resultIndex].setVal(result->getVal(myMatrix[i].getRow(), M.myMatrix[j].getCol()) + myMatrix[i].getVal() * M.myMatrix[j].getVal());
                ++resultIndex;
            }
            
        }
    }
    return result;
}
SparseMatrix* SparseMatrix::Add(SparseMatrix& M){
   // Ensure matrix dimensions are valid for addition
    if (noRows != M.noRows || noCols != M.noCols) {
        throw IncompatibleAdditionException();  // Exception for dimension mismatch
    }

    //Create a counter for the resulting number of sparse values
    int tempNoNonSparseValues = 0;

    //Loop through to check how many non-sparse values there are
    for (int i = 0; i < noRows; i++) {
        for (int j = 0; j < noCols; j++) {
            int sum = nonSparseMatrix[i][j] + M.nonSparseMatrix[i][j];
            //Check if the value is equal to the common value. If not, increment the counter
            if (sum != commonValue) {
                tempNoNonSparseValues++; 
            }
        }
    }

    //Create a new sparse matrix object for the result
    SparseMatrix* result = new SparseMatrix(noRows, noCols, commonValue, tempNoNonSparseValues);

    //Fill the result matrix
    int sparseIndex = 0;
    for (int i = 0; i < noRows; i++) {
        for (int j = 0; j < noCols; j++) {
            int sum = nonSparseMatrix[i][j] + M.nonSparseMatrix[i][j];
            if (sum != commonValue) {
                //Insert the sum into the sparse matrix
                result->myMatrix[sparseIndex].setRow(i);
                result->myMatrix[sparseIndex].setCol(j);
                result->myMatrix[sparseIndex].setVal(sum);
                sparseIndex++;
            }
        }
    }
    return result;
}
ostream& operator<<(ostream& s, const SparseMatrix& sm){
    for(int i = 0; i < sm.getNoNonSparseValues(); ++i){
        s << sm.getMyMatrix()[i];
    }
    return s;
}
void SparseMatrix::displayMatrix(){
    toMatrix();
    //Loop through each row followed by each column to display the matrix
    for(int i = 0; i < noRows; ++i){
        for(int j = 0; j < noCols; ++j){
            cout << nonSparseMatrix[i][j] << " ";
        }
        cout << endl;
    }
}
//getters
int SparseMatrix::getNoRows() const{
    return noRows;
}
int SparseMatrix::getNoCols() const{
    return noCols;
}
int SparseMatrix::getCommonValue() const{
    return commonValue;
}
int SparseMatrix::getNoNonSparseValues() const{
    return noNonSparseValues;
}
SparseRow* SparseMatrix::getMyMatrix() const{
    return myMatrix;
}
int SparseMatrix::getVal(int r, int c) const{
    for(int i = 0; i < noNonSparseValues; ++i){
        //If the row and column match, return the value
        if(myMatrix[i].getRow() == r && myMatrix[i].getCol() == c){
            return myMatrix[i].getVal();
        }
    }
    return commonValue;
}
//setters
void SparseMatrix::setNoRows(int r){
    noRows = r;
}
void SparseMatrix::setNoCols(int c){
    noCols = c;
}
void SparseMatrix::setCommonValue(int cv){
    commonValue = cv;
}
void SparseMatrix::setNoNonSparseValues(int n){
    noNonSparseValues = n;
}
//helper methods
void SparseMatrix::toMatrix(){
    //Create a new matrix to store the sparse matrix in matrix form
    nonSparseMatrix = new int*[noRows];
    for(int i = 0; i < noRows; ++i){
        //Create the dimensions of the matrix
        nonSparseMatrix[i] = new int[noCols];
    }
    for(int i = 0; i < noRows; ++i){
        for(int j = 0; j < noCols; ++j){
            //Fill the matrix with the common value
            nonSparseMatrix[i][j] = commonValue;
        }
    }
    for(int i = 0; i < noNonSparseValues; ++i){
        //Fill the matrix with the sparse values
        nonSparseMatrix[myMatrix[i].getRow()][myMatrix[i].getCol()] = myMatrix[i].getVal();
    }
}
void SparseMatrix::readFile(SparseMatrix& M){
    //Read in the matrix from the file
    for(int i = 0; i < M.noRows; ++i){
        for(int j = 0; j < M.noCols; ++j){
            int temp;
            cin >> temp;
            //If the value is not the common value, add it to the sparse matrix
            if(temp != M.commonValue){
                M.myMatrix[M.matrixIndex].setRow(i);
                M.myMatrix[M.matrixIndex].setCol(j);
                M.myMatrix[M.matrixIndex].setVal(temp);
                ++M.matrixIndex;
            }
        }
    }

}

int main () {
    int n, m, cv, noNSV;
    SparseMatrix* temp;

    // Input for the first matrix
    cin >> n >> m >> cv >> noNSV;
    SparseMatrix* firstOne = new SparseMatrix(n, m, cv, noNSV);

    // Read the first matrix
    firstOne->readFile(*firstOne);

    // Input for the second matrix
    cin >> n >> m >> cv >> noNSV;
    SparseMatrix* secondOne = new SparseMatrix(n, m, cv, noNSV);

    // Read the second matrix
    secondOne->readFile(*secondOne);

    // First matrix in sparse format
    cout << "First one in sparse matrix format" << endl;
    cout << *firstOne;

    // Transpose of the first matrix in sparse format
    cout << "After transpose" << endl;
    temp = firstOne->Transpose();
    cout << *temp;

    // First matrix in matrix format
    cout << "First one in matrix format" << endl;
    firstOne->displayMatrix();

    // Second matrix in sparse format
    cout << "Second one in sparse matrix format" << endl;
    cout << *secondOne;

    // Transpose of the second matrix in sparse format
    cout << "After transpose" << endl;
    temp = secondOne->Transpose();
    cout << *temp;

    // Second matrix in matrix format
    cout << "Second one in matrix format" << endl;
    secondOne->displayMatrix();

    // Matrix addition result
    cout << "Matrix addition result" << endl;
    //Check if the matrices are compatible for addition
    try{
        temp = firstOne->Add(*secondOne);
        temp->displayMatrix();
    }
    catch(IncompatibleAdditionException e){
        cout << e.what() << endl;
    }

    // Matrix multiplication result
    cout << "Matrix multiplication result" << endl;
    //Check if the matrices are compatible for multiplication
    try{
        temp = firstOne->Multiply(*secondOne);
        temp->displayMatrix();
    }
    catch(IncompatibleMultiplicationException e){
        cout << e.what() << endl;
    }

    // Clean up memory
    delete firstOne;
    delete secondOne;
    delete temp;

    return 0;
}

/*

LLM and GitHub Copilot Usage Documentation:

Prompt #1: 
What is the relation between SparseRow and SparseMatrix?
Rationale: 
I was able to clearly understand what the SparseRow class was used for, but was not able to grasp the relation
 to SparseMatrix. I used this prompt to help me understand the relation between the two classes. Copilot was able to explain 
 how the two were likely to be used alongside each other. This made the implementation of the SparseMatrix class much easier, 
 as I was able to understand how the functions within SparseRow would be used within the functions of the SparseMatrix class.

Prompt #2: 
What does it mean to transpose a matrix and how can this be implemented in code?
Rationale: 
After following the dictated layout for the SparseMatrix class, I was unsure what it meant to transpose a 
matrix. I used this prompt to help answer that question. Copilot was able to answer this question very clearly and 
gave me a good understanding of what it means to transpose a matrix. I was able to use this information to implement 
the Transpose() function in the SparseMatrix class.

Prompt #3: 
How does addition function between two sparse matrices?
Rationale: 
While I understood why matrix addition was important, I was not sure of the specifics. I used this prompt to 
understand how two matrices are added together and how to structure a subsequent function. Copilot was able to explain 
matrix addition very clearly, giving me a solid guideline for how to move forward with the implementation of the Add() 
function in the SparseMatrix class.

Prompt #4: 
Implement a condition to check if two matrices are compatible for addition.
Rationale: 
Through my implementation of the Add() function, I realized that I needed to check if the two matrices were 
compatible to be added together in the first place. I was not sure, however, what conditions would need to be met to ensure
compatibility. I used this prompt to help me understand what conditions would need to be met. Copilot was able to explain that, 
for matrix addition, the two matrices must have the same dimensions. I was then easily able to implement the checker for 
compatibility.

Prompt #5: 
Implement a nested for-loop to count the total number of non sparse values when adding two sparse matrices.
Rationale: 
Whilst implementing the Multiply() function, I realized that I needed to count the total number of non-sparse values, 
so that I knew how many elements to allocate for the result matrix. I was not sure how to implement this, so I used this prompt. 
Copilot gave me a clear blueprint to use, which I was then able to alter and improve upon to fit my specific needs. While I still 
had to change the code to fit into the Multiply() function, this prompt gave me a good initial starting point.

Prompt #6: 
How do we account for the common value when multiplying two sparse matrices?
Rationale: 
While attempting to multiply two sparse matrices, I became confused by the common value and whether or not it needed to 
be accounted for. After using this prompt, my view of the Multiply() function became much clearer. Copilot was able to explain that, 
because the matrix was already in sparse form, the common value would not affect the multiplication. This greatly simplified how 
I was viewing the function, and made the implementation much easier.

Prompt #7: 
Define a custom exception class for instances when two matrices are incompatible for multiplication.
Rationale: 
While I knew that I needed to check if two matrices were compatible, I did not know how to create a custom exception class 
in C++. This prompt made it very clear how to do so. Copilot gave a general blueprint for how to create a custom exception class, 
which I then altered to fit my specific needs of multiplication incompatibility. 

Prompt #8: 
Why is my matrix multiplication output slightly off in some cases?
Rationale: 
After implementing the Multiply() function, I noticed that the output was slightly off in a small number of cases. 
After being unable to find the issue, I asked Copilot. It quickly gave me a list of possibilities as to why this was occuring, 
as well as a revised version of the Multiply() function that would likely fix the issue. After reading the suggestions, I 
quickly found the issue and resolved it. 

Prompt #9: 
How should the sparse matrix be displayed?
Rationale: 
While a normal matrix is fairly simple to display, I was unsure of how to do so for a sparse matrix. I used this prompt 
to see what the expected way to display a sparse matrix was. Copilot gave me a clear explanation of how to display a sparse matrix, 
which I was then able to implement into my code. 

Prompt #10: 
Why are there extra newline characters between parts of my output?
Rationale: 
I noticed that, after finalizing my code, there were extra lines between parts of the output. I spent time looking 
for the issue, but was unsuccessful. I used this prompt to ask Copilot why this was happening. Copilot was able to give me a 
list of possibilities, which allowed me to find the issue. I was then able to revise my code and remove the extra newlines, 
which corrected my output to match what we were expected to have. 

Incremental Development: 
Copilot was extremely helpful in the development of my code, specifically as a guide for how to 
begin implementation of various functions. After reading through the assignment description, my understanding of sparse 
matrices was fairly weak. I used Copilot to greatly improve my general understanding of the topic. This was largely done by 
asking questions about the implementation of various aspects of the SparseMatrix class. While Copilot was not used to write any 
functions for me, having a general idea of how to move forward was greatly beneficial to not only my understanding of the functions 
but my efficiency in writing them and the quality of my code. My use of Copilot was mostly for generating a guideline for how to 
move forward, but I also used it on a few occasions for refining my code as well. A few small issues arose throughout, and 
Copilot was very useful in finding and helping eliminate these issues. For the most part, Copilot was used to aid in the 
implementation of the Transpose(), Multiply(), and Add() functions, as well as a few smaller issues. 

----------------------------------------------------------------------------------------------------------------------------

Debugging and Testing Plan: 

Specific Tests:
While implementing various functions in my code, I made sure to test them throughout. This makes debugging much easier, as I 
am able to work through issues as I go rather than attempting to run everything for the first time at the end. If I had done that, 
I would have run into a massive number of issues that would have been very difficult to work through all at once. To test throughout, 
I created separate files of data on top of the ones given to us in the assignment. This allowed me to test matrices with different 
dimensions, common values other than 0, and higher and lower numbers of sparse values. This was very helpful for debugging, especially 
in the case of the Add() and Multiply() functions. These functions were where I had the majority of my issues, and working through 
them while testing likely made my experience much easier than it would have been otherwise. 

Issues and Resolutions:
The majority of my issues arose from the implementation of the Add() and Multiply() functions. These issues were largely due to 
improper use of arrays throughout. Arrays were the primary data structure used in this code, and keeping track of what was happening 
throughout proved to be a challenge. One issue was that I was not checking for the compatibility of matrices before attempting to add
 or multiply them. I fixed this by creating exceptions and using try-catch statements. I also had issued with the output of the
 Multiply() function. I found that my numbers were slightly off in different cases. I discovered that this was caused by not accounting 
for some cases. Using Copilot to help explain as well as several new inputs, I was able to resolve the issue. The final large issue 
I ran into was adding matrices in sparse form. I found that it was much simpler to do addition in normal matrix form, and then convert 
to sparse form. I solved this by creating a few helper methods, namely toMatrix() and readFile(). These methods allowed me to easily 
convert between the two forms, and made the implementation of the Add() function much easier.

Verification: 
My primary mode of verification was taking inputs, both from the assignment and from my own files, and comparing them to the expected 
outputs. I began with the five input examples on the assignment, and compared them to the expected output. I then used my own test cases, 
making sure to use different common values and dimensions. This ensured that my code was able to handle common values other than 0 as well 
as matrices with different dimensions or number of sparse values. After testing several different sets of numbers, I was satisfied that my 
code would be able to handle not only the input from the assignment, but any other inputs as well. 

*/