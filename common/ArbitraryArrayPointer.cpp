//////////////////////////////////////////////////////////////////////////
// A class that can point to any Array<T>, to simplify class Particles<d>
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ArbitraryArrayPointer.h"

void ArbitraryArrayPointer::Resize(const int size)
{
	data_ptr->Resize(size);
}

void ArbitraryArrayPointer::Reserve(const int size)
{
	data_ptr->Reserve(size);
}

int ArbitraryArrayPointer::Add_Element(void)
{
	return data_ptr->Add_Element();
}

int ArbitraryArrayPointer::Add_Elements(const int n)
{
	return data_ptr->Add_Elements(n);
}

void ArbitraryArrayPointer::Copy_Element(const int idx, const ArbitraryArrayPointer& src_array, const int src_idx)
{
	data_ptr->Copy_Element(idx, src_array.data_ptr, src_idx);
}

int ArbitraryArrayPointer::Delete_Elements(const Array<bool>& is_deleted)
{
	return data_ptr->Delete_Elements(is_deleted);
}

bool ArbitraryArrayPointer::Write_Binary(const std::string& file_name)
{
	return data_ptr->Write_Binary(file_name);
}

bool ArbitraryArrayPointer::Read_Binary(const std::string& file_name)
{
	return data_ptr->Read_Binary(file_name);
}

