LOCAL_PATH:= $(call my-dir)

ROOT_DIR:=../../../
SRC_DIR:=$(ROOT_DIR)src/
INCLUDE_DIR:=jni/../../../include/

include $(CLEAR_VARS)
LOCAL_C_INCLUDES += $(SRC_DIR) $(INCLUDE_DIR)
LOCAL_MODULE    := ezsift
LOCAL_CFLAGS := -O3 -ffast-math -fPIE
LOCAL_SRC_FILES := $(SRC_DIR)/image_utility.cpp \
                   $(SRC_DIR)/ezsift.cpp
include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)
LOCAL_SRC_FILES:=$(ROOT_DIR)/examples/feature_extract/feature_extract.cpp 
LOCAL_C_INCLUDES += $(SRC_DIR) $(INCLUDE_DIR)
LOCAL_MODULE:= feature_extract
LOCAL_CFLAGS := -O3 -ffast-math -fPIE
LOCAL_STATIC_LIBRARIES := ezsift
include $(BUILD_EXECUTABLE)

include $(CLEAR_VARS)
LOCAL_SRC_FILES := $(ROOT_DIR)/examples/image_match/image_match.cpp 
LOCAL_C_INCLUDES += $(SRC_DIR) $(INCLUDE_DIR)
LOCAL_MODULE:= image_match
LOCAL_STATIC_LIBRARIES := ezsift
LOCAL_CFLAGS := -O3 -ffast-math -fPIE
LOCAL_LDFLAGS := -fPIE -pie 
include $(BUILD_EXECUTABLE)