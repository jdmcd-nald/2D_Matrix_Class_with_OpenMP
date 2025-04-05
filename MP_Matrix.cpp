#include <iostream>
#include <fstream>
#include <filesystem>
#include <mdspan>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <utility>
#include <string>

#include <array>
#include <iostream>
#include <any>

#include <cassert>

//#include <omp.h>

template <typename T, int N, int M>
int length(T(&)[N][M])
{
    return N;
}

class Matrix2D
{
private:
    // You might be tempted to make m_arr a reference to an ArrayFlat2d,
    // but this makes the view non-copy-assignable since references can't be reseated.
    // Using std::reference_wrapper gives us reference semantics and copy assignability.
    std::vector<long double> my2DArray; // 3x3 2D array
    int m_row;
    int m_column;

public:
    Matrix2D(int row,int col,std::vector<long double> vec) : m_row(row),m_column(col),my2DArray(vec) {}

    // Get element via single subscript (using operator[])
    //T& operator[](int i) { return m_arr[i]; }
    //const T& operator[](int i) const { return m_arr[i]; }

    // Get element via 2d subscript (using operator(), since operator[] doesn't support multiple dimensions prior to C++23)
    //T& operator[](int row, int col) { return m_arr[static_cast<std::size_t>(row * cols() + col)]; }
    //const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }

    
    
    // in C++23, you can uncomment these since multidimensional operator[] is supported
//    T& operator[](int row, int col) { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
//    const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
    
    void replace_data(Matrix2D& arr)
    {
        assert((*this).rows() == arr.rows());
        assert(       (*this).cols() == arr.cols());
        
        for(int k = 0;k<(*this).rows()*(*this).cols();k++)
        {
            ((*this).data())[k] = ((arr).data())[k];
        }
    }
    
    void scale(long double scaler)
    {
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_2 = scaler *((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2;
            }
        }
    }
    
    void addition(Matrix2D& arr)
    {
        assert((*this).cols() == arr.cols());
        assert(       (*this).rows() == arr.rows());
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_1+test_2;
            }
        }
    }
    
    void subtraction(Matrix2D& arr)
    {
        assert((*this).rows() == arr.rows());
        assert(       (*this).cols() == arr.cols());
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2-test_1;
            }
        }
    }
    
    Matrix2D multiply(Matrix2D& arr)
    {
        
        assert((*this).cols() == arr.rows());
        
        std::vector<long double> output_vector;
        output_vector.reserve((*this).rows()*arr.cols());
        long double holder = 0;
        for(int j = 0;j<(*this).rows();j++)
        {
            for(int k = 0;k<(arr).cols();k++)
            {
                holder = 0;
                for(int i = 0;i<(*this).cols();i++)
                {
                    long double left = ((*this).data())[j*((*this).cols())+i];
                    long double right = (arr.data())[i*(arr.cols())+k];
                    holder += left*right;
                }
                output_vector.push_back(holder);
            }
        }
        Matrix2D output((*this).rows(),arr.cols(),output_vector);
        return output;
    }
    
    Matrix2D output_identity(int size)
    {
        std::vector<long double> output_vector;
        output_vector.reserve(size*size);
        for(int j = 0;j<size;j++)
        {
            for(int k = 0;k<size;k++)
            {
                if(j==k)
                {
                    output_vector.push_back(1);
                }
                else
                {
                    output_vector.push_back(0);
                }
            }
        }
        Matrix2D output(size,size,output_vector);
        return output;
    }
    
    Matrix2D transpose()
    {
        std::vector<long double> output_vector;
        output_vector.reserve((*this).rows()*(*this).cols());
        for(int j = 0;j<(*this).cols();j++)
        {
            for(int k = 0;k<(*this).rows();k++)
            {
                output_vector.push_back(((*this).data())[k*((*this).cols())+j]);
            }
        }
        Matrix2D output((*this).cols(),(*this).rows(),output_vector);
        return output;
    }
    
    
    
    Matrix2D cofactor()
    {
        assert((*this).rows() == (*this).cols());
        
        //std::vector<Matrix2D> cofactor_segment;
        //cofactor_segment.reserve((*this).rows() * (*this).cols());
        std::vector<long double> cofactor_matrix_vector;
        
        if(((*this).rows()) == 2)
        {
            cofactor_matrix_vector.push_back(((*this).data())[3]);
            cofactor_matrix_vector.push_back(-((*this).data())[1]);
            cofactor_matrix_vector.push_back(-((*this).data())[2]);
            cofactor_matrix_vector.push_back(((*this).data())[0]);
            Matrix2D cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
            return cofactor_matrix;
        }
        
        for(int j = 0;j<(*this).cols();j++)
        {
            for(int k = 0;k<(*this).rows();k++)
            {
                std::vector<long double> cofactor_segment_vector;
                cofactor_segment_vector.reserve((((*this).rows())-1) * (((*this).cols())-1));
                
                for(int i = 0;i<(*this).rows();i++)
                {
                    for(int o = 0;o<(*this).cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            cofactor_segment_vector.push_back(((*this).data())[i*((*this).cols())+o]);
                        }
                    }
                }
                
                int u = -1;
                if(k%2 == j%2)
                {
                    u=1;
                }
                
                Matrix2D temp(((*this).rows())-1,((*this).cols())-1,cofactor_segment_vector);
                
                cofactor_matrix_vector.push_back(u*det(temp));
                
                //cofactor_segment.push_back(temp);
            }
        }
        
        Matrix2D cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
        return cofactor_matrix;
    }
    
    double long det(Matrix2D arr)
    {
        std::vector<long double> determinate;
        if((arr).rows()==2){
            return(det_2D(arr));
        }
        
        for(int j = 0;j<arr.cols();j++)
        {
            int k =0;
            {
                std::vector<long double> det_segment_vector;
                det_segment_vector.reserve((((arr).rows())-1) * (((arr).cols())-1));
                
                for(int i = 0;i<arr.rows();i++)
                {
                    for(int o = 0;o<arr.cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            det_segment_vector.push_back(((arr).data())[i*((arr).cols())+o]);
                        }
                    }
                }
                
                int u = -1;
                if(k%2 == j%2){ u=1; }
                
                Matrix2D temp((((arr).rows())-1) ,(((arr).cols())-1),det_segment_vector);
                
                determinate.push_back(u*(((arr).data())[k*((arr).cols())+j])*det(temp));
                
                //cofactor_segment.push_back(temp);
            }
        }
        long double sum = 0;
        for(int i =0;i<determinate.size();i++)
        {
            sum += determinate[i];
        }
        return sum;
    }
    
    Matrix2D inverse(Matrix2D& arr)
    {
        Matrix2D output = ((arr.cofactor()).transpose());
        long double test = det(arr);
        output.scale(1.0/det(arr));
        return output;
    }
    
    long double det_2D(Matrix2D& arr)
    {
        assert(arr.rows() == arr.cols());
        return (arr.data())[0]*(arr.data())[3]-(arr.data())[1]*(arr.data())[2];
    }
    
    void print()
    {
        for( int i =0 ; i<this->rows() ; i++)
        {
            for( int j =0 ; j<this->cols() ; j++)
            {
                std::cout << (this->data())[i*(this->cols())+j]<<"\t";
            }
            std::cout << std::endl;
        }
        std::cout <<this->rows()<<","<<this->cols()<< std::endl;
    }
    
    int rows() const { return m_row; }
    int cols() const { return m_column; }
    std::vector<long double>& data() { return my2DArray; }
};


class Matrix2D_MP
{
private:
    // You might be tempted to make m_arr a reference to an ArrayFlat2d,
    // but this makes the view non-copy-assignable since references can't be reseated.
    // Using std::reference_wrapper gives us reference semantics and copy assignability.
    std::vector<long double> my2DArray; // 3x3 2D array
    int m_row;
    int m_column;

public:
    Matrix2D_MP(int row,int col, std::vector<long double> vec ) : m_row(row),m_column(col),my2DArray(vec) {}

    // Get element via single subscript (using operator[])
    //T& operator[](int i) { return m_arr[i]; }
    //const T& operator[](int i) const { return m_arr[i]; }

    // Get element via 2d subscript (using operator(), since operator[] doesn't support multiple dimensions prior to C++23)
    //T& operator[](int row, int col) { return m_arr[static_cast<std::size_t>(row * cols() + col)]; }
    //const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }

    
    
    // in C++23, you can uncomment these since multidimensional operator[] is supported
//    T& operator[](int row, int col) { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
//    const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
    
    void replace_data(Matrix2D_MP& arr)
    {
        assert((*this).rows() == arr.rows());
        assert((*this).cols() == arr.cols());
        
        int k;
        #pragma omp paralllel for
        for(k = 0;k<(*this).rows()*(*this).cols();k++)
        {
            ((*this).data())[k] = ((arr).data())[k];
        }
    }
    
    void scale(long double scaler)
    {
        //std::vector<long double> output_vector;
        int k,j;
        #pragma omp paralllel for private(j)
        for(k = 0;k<(*this).rows();k++)
        {
            for(j = 0;j<(*this).cols();j++)
            {
                long double test_2 = scaler *((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2;
            }
        }
    }
    
    void addition(Matrix2D_MP& arr)
    {
        assert((*this).cols() == arr.cols());
        assert((*this).rows() == arr.rows());
        //std::vector<long double> output_vector;
        int k,j;
        #pragma omp paralllel for private(j)
        for(k = 0;k<(*this).rows();k++)
        {
            for(j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_1+test_2;
            }
        }
    }
    
    void subtraction(Matrix2D_MP& arr)
    {
        assert((*this).rows() == arr.rows());
        assert((*this).cols() == arr.cols());
        
        //std::vector<long double> output_vector;
        
        int k,j;
        #pragma omp paralllel for private(j)
        for(k = 0;k<(*this).rows();k++)
        {
            for(j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2-test_1;
            }
        }
    }
    
    Matrix2D_MP multiply(Matrix2D_MP& arr)
    {
        assert((*this).cols() == arr.rows());
        
        std::vector<long double> output_vector((*this).rows()*arr.cols(),0);
        //output_vector.reserve((*this).rows()*arr.cols());
        long double holder = 0;
        int j,k,i;
        #pragma omp parallel for private(k,j,i,holder) collapse(2)
        for(j = 0;j<(*this).rows();j++)
        {
            for(k = 0;k<(arr).cols();k++)
            {
                holder = 0;
                for(i = 0;i<(*this).cols();i++)
                {
                    long double left = ((*this).data())[j*((*this).cols())+i];
                    long double right = (arr.data())[i*(arr.cols())+k];
                    holder += left*right;
                }
                
                //output_vector.push_back(holder);
                #pragma omp atomic
                output_vector[j*(arr).cols()+k]=holder;
            }
        }
        Matrix2D_MP output((*this).rows(),arr.cols(),output_vector);
        return output;
    }
    
    Matrix2D_MP output_identity(int size)
    {
        std::vector<long double> output_vector(size*size,0);
        //output_vector.reserve(size*size);
        for(int j = 0;j<size;j++)
        {
            //if(j==k)
            //{
            output_vector[j*size+j]=1;
            //}
            //else
            //{
            //    output_vector.push_back(0);
            //}
        }
        Matrix2D_MP output(size,size,output_vector);
        return output;
    }
    
    Matrix2D_MP transpose()
    {
        std::vector<long double> output_vector(((*this).rows()*(*this).cols()),0);
        //output_vector.reserve((*this).rows()*(*this).cols());
        int j=0;
        int k=0;
        #pragma omp parallel for private(j,k) collapse(2)
        for(j = 0;j<(*this).cols();j++)
        {
            for(k = 0;k<(*this).rows();k++)
            {
                output_vector[j*(*this).rows()+k]=(((*this).data())[k*((*this).cols())+j]);
            }
        }
        Matrix2D_MP output((*this).cols(),(*this).rows(),output_vector);
        return output;
    }
    
    
    
    Matrix2D_MP cofactor()
    {
        assert((*this).rows() == (*this).cols());
        
        //std::vector<Matrix2D_MP> cofactor_segment;
        //cofactor_segment.reserve((*this).rows() * (*this).cols());
        std::vector<long double> cofactor_matrix_vector;
        
        if(((*this).rows()) == 2)
        {
            cofactor_matrix_vector.push_back(((*this).data())[3]);
            cofactor_matrix_vector.push_back(-((*this).data())[1]);
            cofactor_matrix_vector.push_back(-((*this).data())[2]);
            cofactor_matrix_vector.push_back(((*this).data())[0]);
            Matrix2D_MP cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
            return cofactor_matrix;
        }
        
        int j;
        int k;
        std::vector<long double> cofactor_segment_vector(((((*this).rows())-1) * (((*this).cols())-1)),0);
        #pragma omp parallel for private(j,k) public(cofactor_segment_vector) collapse(2)
        for(j = 0;j<(*this).cols();j++)
        {
            for(k = 0;k<(*this).rows();k++)
            {
                for(int i = 0;i<(*this).rows();i++)
                {
                    for(int o = 0;o<(*this).cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            cofactor_segment_vector[i*(this->cols())+o]=(((*this).data())[i*((*this).cols())+o]);
                        }
                    }
                }
                int u = -1;
                if(k%2 == j%2)
                {
                    u=1;
                }
                Matrix2D_MP temp(((*this).rows())-1,((*this).cols())-1,cofactor_segment_vector);
                long double temporary = (u*det(temp));
                #pragma omp atomic
                cofactor_matrix_vector[j*((*this).rows())+k]=temporary;
                //cofactor_segment.push_back(temp);
            }
        }
        
        Matrix2D_MP cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
        return cofactor_matrix;
    }
    
    double long det(Matrix2D_MP arr)
    {
        
        if((arr).rows()==2){
            return(det_2D(arr));
        }
        
        std::vector<long double> det_segment_vector(((((arr).rows())-1) * (((arr).cols())-1)),0);
        //det_segment_vector.reserve((((arr).rows())-1) * (((arr).cols())-1));
        std::vector<long double> determinate(arr.cols(),0);
        for(int j = 0;j<arr.cols();j++)
        {
            int k =0;
            {
                //std::vector<long double> det_segment_vector;
                //det_segment_vector.reserve((((arr).rows())-1) * (((arr).cols())-1));
                
                for(int i = 0;i<arr.rows();i++)
                {
                    for(int o = 0;o<arr.cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            det_segment_vector[i*arr.cols()+o]=(((arr).data())[i*((arr).cols())+o]);
                        }
                    }
                }
                
                int u = -1;
                if(k%2 == j%2){ u=1; }
                
                Matrix2D_MP temp((((arr).rows())-1) ,(((arr).cols())-1),det_segment_vector);
                
                determinate[j]=(u*(((arr).data())[k*((arr).cols())+j])*det(temp));
                
                //cofactor_segment.push_back(temp);
            }
        }
        long double sum = 0;
        for(int i =0;i<determinate.size();i++)
        {
            sum += determinate[i];
        }
        return sum;
    }
    
    Matrix2D_MP inverse(Matrix2D_MP& arr)
    {
        Matrix2D_MP output = ((arr.cofactor()).transpose());
        long double test = det(arr);
        output.scale(1.0/det(arr));
        return output;
    }
    
    long double det_2D(Matrix2D_MP& arr)
    {
        assert(arr.rows() == arr.cols());
        return (arr.data())[0]*(arr.data())[3]-(arr.data())[1]*(arr.data())[2];
    }
    
    void print()
    {
        for( int i =0 ; i<this->rows() ; i++)
        {
            for( int j =0 ; j<this->cols() ; j++)
            {
                std::cout << (this->data())[i*(this->cols())+j]<<"\t";
            }
            std::cout << std::endl;
        }
        std::cout <<this->rows()<<","<<this->cols()<< std::endl;
    }
    
    int rows() const { return m_row; }
    int cols() const { return m_column; }
    std::vector<long double>& data() { return my2DArray; }
};

int main(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";
    
    Matrix2D cofactor_test_2(2,2,{1,2,
                                3,4});
    Matrix2D_MP cofactor_test_22(2,2,
                                 {1,2,
                                 3,4});
    
    cofactor_test_2.print();
    cofactor_test_22.print();
    
    cofactor_test_2.cofactor().print();
    cofactor_test_2.cofactor().print();
    
    Matrix2D cofactor_test_3(3,3,
        {1,2,3,
         4,5,6,
         7,8,9});
    
    
    
    Matrix2D cofactor_test_33(3,3,
        {1,2,3,
         2,5,4,
         3,3,2});
    
    cofactor_test_3.cofactor().print();
    cofactor_test_33.cofactor().print();
    
    
    Matrix2D cofactor_test_4(4,4,
        {1,2,3,4,
         2,4,1,1,
         3,4,1,2,
         4,2,1,3});
    
    Matrix2D cofactor_test_44(4,4,
        {1,2,3,4,
         2,4,1,1,
         3,4,1,2,
         4,2,1,3});
    
    cofactor_test_4.cofactor().print();
    cofactor_test_44.cofactor().print();
    
    Matrix2D cofactor_test_5(4,4,
        {9,9,3,4,
         9,8,5,19,
         4,5,7,2,
         4,2,9,3});
    
    Matrix2D cofactor_test_55(4,4,
        {9,9,3,4,
         9,8,5,19,
         4,5,7,2,
         4,2,9,3});
    
    
    cofactor_test_5.inverse(cofactor_test_5).print();
    cofactor_test_55.inverse(cofactor_test_55).print();
    
    Matrix2D test(3,3,
        {6,2,4,
        -1,4,3,
        -2,9,3});
    
    Matrix2D test1(3,1,
        {4,
        -2,
         1});
    
    Matrix2D testt(3,3,
        {2,3,2,
         2,3,1,
         2,1,0});
    
    Matrix2D testt1(3,2,
        {3,1,
         3,2,
         1,0});
    
    test.multiply(test1);
    testt.multiply(testt1);
    auto id = testt.output_identity(3);
    
    Matrix2D test_mp(3,3,
        {6,2,4,
        -1,4,3,
        -2,9,3});
    
    Matrix2D test1_mp(3,1,
        {4,
        -2,
         1});
    
    Matrix2D testt_mp(3,3,
        {2,3,2,
         2,3,1,
         2,1,0});
    
    Matrix2D testt1_mp(3,2,
        {3,1,
         3,2,
         1,0});
    
    test_mp.multiply(test1_mp);
    testt_mp.multiply(testt1_mp);
    auto id_mp = testt_mp.output_identity(3);
    
    test.print();
    test_mp.print();
    testt.print();
    testt_mp.print();
    id.print();
    id_mp.print();
    
    
    
    Matrix2D transpose_test(2,3,
        {1,2,3,
         4,5,6});
    
    (transpose_test.transpose()).print();
    
    Matrix2D_MP transpose_test_mp(2,3,
        {1,2,3,
         4,5,6});
    
    (transpose_test_mp.transpose()).print();
    
    Matrix2D transpose_test_2(3,3,
        {1,2,3,
        4,5,6,
        7,8,9});
    
    (transpose_test_2.transpose()).print();
    
    Matrix2D_MP transpose_test_2_mp(3,3,
        {1,2,3,
        4,5,6,
        7,8,9});
    
    (transpose_test_2_mp.transpose()).print();
    
    
    Matrix2D addition_test(2,2,
        {8,5,
        2,3});
    Matrix2D addition_test_2(2,2,
        {9,5,
        1,2});
    addition_test.addition(addition_test_2);
    addition_test.print();
    
    Matrix2D addition_test_mp(2,2,
        {8,5,
        2,3});
    Matrix2D addition_test_2_mp(2,2,
        {9,5,
        1,2});
    addition_test_mp.addition(addition_test_2_mp);
    addition_test_mp.print();
   
   
    /*
    Matrix2D cofactor_test(3,3,
        {11,12,13,
         21,22,23,
         31,32,33});
    */
    return 0;
}
