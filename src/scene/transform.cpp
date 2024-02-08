#include "transform.h"
#include <iostream>

Mat4 Transform::local_to_parent() const {
	return Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
}

Mat4 Transform::parent_to_local() const {
	return Mat4::scale(1.0f / scale) * rotation.inverse().to_mat() * Mat4::translate(-translation);
}

Mat4 Transform::local_to_world() const {
	// A1T1: local_to_world
	//don't use Mat4::inverse() in your code.
    Mat4 local_to_world_transformation = Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
    if (std::shared_ptr< Transform > spt = parent.lock()) {
        //case where transform has a parent
        Mat4 parent_to_world_transformation = spt->local_to_world();
        return parent_to_world_transformation * local_to_world_transformation;
    }
    else {
        //case where transform doesn't have a parent
        return local_to_world_transformation;
    }
//    return Mat4::I; //<-- wrong, but here so code will compile
}

Mat4 Transform::world_to_local() const {
	// A1T1: world_to_local
	//don't use Mat4::inverse() in your code.
    Mat4 world_to_local_transformation = Mat4::scale(1.0f/scale) * Mat4::transpose(rotation.to_mat()) * Mat4::translate(-translation);
    if (std::shared_ptr< Transform > spt = parent.lock()) {
        //case where transform has a parent
        Mat4 parent_to_local_transformation = spt->world_to_local();
        return world_to_local_transformation * parent_to_local_transformation;
    }
    else {
        //case where transform doesn't have a parent
        return world_to_local_transformation;
    }
//	return Mat4::I; //<-- wrong, but here so code will compile
}

bool operator!=(const Transform& a, const Transform& b) {
	return a.parent.lock() != b.parent.lock() || a.translation != b.translation ||
	       a.rotation != b.rotation || a.scale != b.scale;
}
