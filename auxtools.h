#ifndef __auxtools_h__
#define __auxtools_h__


template<typename T, size_t SIZE>
void dump_table_to_file(simpleComplex<T> *tbl, const char *filename);

template<typename T, size_t SIZE>
void dump_table_to_file(simpleComplex<T> *tbl, const char *matname, const char *filename);


template<typename T>
void dump_uCSRMatrix_to_file(const uCSRMatrix< simpleComplex<T> > &tmat, const char *matname, const char *filename);


#endif // __auxtools_h__