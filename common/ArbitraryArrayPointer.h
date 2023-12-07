//////////////////////////////////////////////////////////////////////////
// A class that can point to any Array<T>, to simplify class Particles<d>
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArbitraryArrayPointer_h__
#define __ArbitraryArrayPointer_h__
#include <iostream>
#include "Common.h"
#include "DataIO.h"

class ArrayPointerContainerBase {
public:
	virtual void Resize(const int size) = 0;
	virtual void Reserve(const int size) = 0;
	virtual int Add_Element(void) = 0;
	virtual int Add_Elements(const int n) = 0;
	virtual void Copy_Element(const int idx, const std::shared_ptr<ArrayPointerContainerBase> src, const int src_idx) = 0;
	virtual int Delete_Elements(const Array<bool>& is_deleted) = 0;//return size after delete
	virtual bool Read_Binary(const std::string& file_name) = 0;
	virtual bool Write_Binary(const std::string& file_name) = 0;
};

template<class T>
class ArrayPointerContainerDerived :public ArrayPointerContainerBase {
public:
	ArrayPtr<T> array_ptr;
	virtual void Resize(const int size) { if (size == 0)array_ptr->clear(); else array_ptr->resize(size); }
	virtual void Reserve(const int size) { array_ptr->reserve(size); }
	virtual int Add_Element(void) { array_ptr->push_back(Zero<T>()); return (int)array_ptr->size() - 1; }
	virtual int Add_Elements(const int n) { array_ptr->resize(array_ptr->size() + n, Zero<T>()); return (int)array_ptr->size() - n; }
	virtual void Copy_Element(const int idx, const std::shared_ptr<ArrayPointerContainerBase> src, const int src_idx) {
		auto derived_src = std::dynamic_pointer_cast<ArrayPointerContainerDerived<T>>(src);
		if (derived_src == nullptr) {
			std::cerr << "ArrayPointerContainerDerived::Copy_Element Error: invalid type\n";
			exit(0);
		}
		const ArrayPtr<T> derived_ptr = derived_src->array_ptr;
		(*array_ptr)[idx] = (*derived_ptr)[src_idx];
	}
	virtual int Delete_Elements(const Array<bool>& is_deleted) {
		int deleted_size = 0;
		for (int i = 0; i < array_ptr->size(); i++) {
			if (is_deleted[i]) continue;
			(*array_ptr)[deleted_size] = (*array_ptr)[i];
			deleted_size++;
		}
		array_ptr->resize(deleted_size);
		return deleted_size;
	}
	virtual bool Read_Binary(const std::string& file_name) { return BinaryDataIO::Read_Array(file_name, *array_ptr); }
	virtual bool Write_Binary(const std::string& file_name) { return BinaryDataIO::Write_Array(file_name, *array_ptr); }
	ArrayPointerContainerDerived(ArrayPtr<T> _ptr) :array_ptr(_ptr) { }
};


class ArbitraryArrayPointer {
public:
	std::shared_ptr<ArrayPointerContainerBase> data_ptr;
	void Resize(const int size);
	void Reserve(const int size);
	int Add_Element(void);
	int Add_Elements(const int n);
	void Copy_Element(const int idx, const ArbitraryArrayPointer& src_array, const int src_idx);
	int Delete_Elements(const Array<bool>& is_deleted);//return size after delete
	bool Write_Binary(const std::string& file_name);
	bool Read_Binary(const std::string& file_name);
	ArbitraryArrayPointer() :data_ptr(nullptr) {}
	template<class T>
	ArbitraryArrayPointer(ArrayPtr<T> arr_ptr) {
		data_ptr = std::make_shared<ArrayPointerContainerDerived<T>>(arr_ptr);
	}
};

#endif
