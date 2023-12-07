#pragma once

#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

namespace YAML {
#define Partially_Specialize_Eigen_Type(Type) \
	template <> \
	struct convert<Eigen::##Type> \
	{ \
		static Node encode(const Eigen::##Type &rhs) \
		{ \
			Node node; \
			for (int i = 0; i < rhs.size(); i++) { \
				node.push_back(rhs(i)); \
			} \
			return node; \
		} \
		static bool decode(const Node &node, Eigen::##Type &rhs) \
		{ \
			if (!node.IsSequence() || node.size() != rhs.size()) { \
				return false; \
			} \
			for (int i = 0; i < rhs.size(); i++) { \
				rhs(i) = node[i].as<Eigen::##Type::Scalar>(); \
			} \
			return true; \
		} \
	}
	Partially_Specialize_Eigen_Type(Vector2i);
	Partially_Specialize_Eigen_Type(Vector3i);
	Partially_Specialize_Eigen_Type(Vector2d);
	Partially_Specialize_Eigen_Type(Vector3d);
#undef Partially_Specialize_Eigen_Type
}
