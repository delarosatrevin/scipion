/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
#include "xmipp_MetaData.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/metadata.h>
#include <data/metadata_label.h>
#include <data/metadata_extension.h>
#include <classification/analyze_cluster.h>
#include "xmipp_InternalData.h"

int debug = 0;
/*
JNIEXPORT void JNICALL Java_xmipp_MetaData_storeIds
(JNIEnv *env, jclass cls) {
 MetaData_peerId = env->GetFieldID(cls, "peer", "J");
}
*/
JNIEXPORT void JNICALL Java_xmipp_MetaData_create
(JNIEnv *env, jobject jobj)
{
    MetaData * md = new MetaData();
    //env->SetLongField(jobj, MetaData_peerId, (long)md);
    STORE_PEER_ID(jobj, md);
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_destroy
(JNIEnv *env, jobject jobj)
{
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    delete md;
    md = NULL;
    //env->SetLongField(jobj, MetaData_peerId, (long)md);
    STORE_PEER_ID(jobj, md);
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_read_1
(JNIEnv *env, jobject jobj, jstring filename)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        const char * fn = env->GetStringUTFChars(filename, false);

        md->read(fn);

        env->ReleaseStringUTFChars(filename, fn);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jint JNICALL Java_xmipp_MetaData_size(JNIEnv *env, jobject jobj)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        return md->size();
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return 0;
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_setColumnFormat(JNIEnv *env, jobject jobj, jboolean format)
{
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    std::string msg = "";

    try
    {
        md->setColumnFormat(format);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_write
(JNIEnv *env, jobject jobj, jstring filename)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        const char * fn = env->GetStringUTFChars(filename, false);

        md->write(fn);

        env->ReleaseStringUTFChars(filename, fn);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_writeBlock
(JNIEnv *env, jobject jobj, jstring filename)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        const char * fn = env->GetStringUTFChars(filename, false);

        md->write(fn, MD_APPEND);

        env->ReleaseStringUTFChars(filename, fn);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_print
(JNIEnv *env, jobject jobj)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        md->write(std::cout);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_containsLabel(JNIEnv *env,
        jobject jobj, jint label)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        return md->containsLabel((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return false;
}

JNIEXPORT jstring JNICALL Java_xmipp_MetaData_label2Str(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";

    try
    {
        String str = MDL::label2Str((MDLabel) label);
        return env->NewStringUTF(str.data());
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jobjectArray JNICALL Java_xmipp_MetaData_getBlocksInMetaDataFile(
    JNIEnv *env, jclass class_, jstring filename)
{
    std::string msg = "";

    try
    {
        std::vector < std::string > blocks;

        const char * fn = env->GetStringUTFChars(filename, false);

        getBlocksInMetaDataFile(fn, blocks);

        env->ReleaseStringUTFChars(filename, fn);

        // Sets array value
        jstring str;
        jobjectArray array = env->NewObjectArray(blocks.size(), env->FindClass(
                                 "java/lang/String"), NULL);

        for (int i = 0; i < blocks.size(); i++)
        {
            str = env->NewStringUTF(blocks[i].data());
            env->SetObjectArrayElement(array, i, str);
            env->DeleteLocalRef(str);
        }

        return array;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jintArray JNICALL Java_xmipp_MetaData_getActiveLabels(JNIEnv *env,
        jobject jobj)
{
    std::string msg = "";
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    try
    {
        std::vector < MDLabel > labels = md->getActiveLabels();

        // Copies vector into array.
        size_t size = labels.size();
        jint *body = new jint[size];
        for (int i = 0; i < size; i++)
        {
            body[i] = labels[i];
        }

        // Sets array value
        jintArray array = env->NewIntArray(size);
        env->SetIntArrayRegion(array, 0, size, body);

        return array;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jclass JNICALL Java_xmipp_MetaData_getLabelType(JNIEnv *env,
        jclass jclass_, jint label)
{
    std::string msg = "";
    jclass class_ = NULL;

    try
    {
        switch (MDL::labelType((MDLabel) label))
        {
        case LABEL_BOOL:
            class_ = env->FindClass("java/lang/Boolean");
            break;
        case LABEL_INT:
            class_ = env->FindClass("java/lang/Integer");
            break;
        case LABEL_LONG:
            class_ = env->FindClass("java/lang/Long");
            break;
        case LABEL_DOUBLE:
            class_ = env->FindClass("java/lang/Double");
            break;
        case LABEL_VECTOR:
        case LABEL_VECTOR_LONG:
        case LABEL_STRING:
            class_ = env->FindClass("java/lang/String");
            break;
        }
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return class_;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isTextFile(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isTextFile((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isMetadata(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isMetadata((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isCtfParam(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isCtfParam((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isImage(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isImage((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isStack(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isStack((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isMicrograph(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isMicrograph((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_isPSD(JNIEnv *env,
        jclass class_, jint label)
{
    std::string msg = "";
    bool result = false;

    try
    {
        result = MDL::isPSD((MDLabel) label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return result;
}

JNIEXPORT jint JNICALL Java_xmipp_MetaData_getValueInt(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        int value;

        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return 0;
}

JNIEXPORT jlong JNICALL Java_xmipp_MetaData_getValueLong
(JNIEnv *env, jobject jobj, jint label, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        size_t value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return 0;
}

JNIEXPORT jdouble JNICALL Java_xmipp_MetaData_getValueDouble(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        double value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return 0;
}

JNIEXPORT jstring JNICALL Java_xmipp_MetaData_getValueString(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    String msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        MDObject obj((MDLabel) label);
        if (md->getValue(obj, objId))
        {
            String str = obj.toString();
            return env->NewStringUTF(str.data());
        }
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_getValueBoolean(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        bool value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return 0;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_setValueInt(JNIEnv *env,
        jobject jobj, jint label, jint value, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        return md->setValue((MDLabel) label, value, objId);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_setValueDouble(JNIEnv *env,
        jobject jobj, jint label, jdouble value, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        return md->setValue((MDLabel) label, value, objId);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return false;
}
JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_setValueString(JNIEnv *env,
        jobject jobj, jint label, jstring value, jlong objId)
{
    std::string msg = "";

    try
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

//        std::cout <<  (MDLabel) label << ": " << MDL_IMAGE << " -> " << MDL::label2Str(MDL_IMAGE) << std::endl;

        const char * strValue = env->GetStringUTFChars(value, false);
        bool result = md->setValueFromStr((MDLabel) label, strValue, objId);

        env->ReleaseStringUTFChars(value, strValue);

        return result;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_MetaData_setValueBoolean(JNIEnv *env,
        jobject jobj, jint label, jboolean value, jlong objId)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        return md->setValue((MDLabel) label, (bool) value, objId);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return false;
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_MetaData_getStatistics(JNIEnv *env,
        jobject jobj, jboolean applyGeo)
{
    std::string msg = "";
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    try
    {
        double avg, stddev, min, max;
        getStatistics((*md), avg, stddev, min, max, applyGeo);

        // Copies vector into array.
        double statistics[4];
        statistics[0] = min;
        statistics[1] = max;
        statistics[2] = avg;
        statistics[3] = stddev;

        // Sets array value
        jdoubleArray array = env->NewDoubleArray(4);
        env->SetDoubleArrayRegion(array, 0, 4, statistics);

        return array;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_MetaData_getColumnValues(JNIEnv *env,
        jobject jobj, jint label)
{
    std::string msg = "";
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    try
    {
        std::vector<double> values;
        md->getColumnValues((MDLabel) label, values);

        // Copies vector into array.
        size_t size = values.size();
        jdouble *body = new jdouble[size];
        for (int i = 0; i < size; i++)
        {
            body[i] = values[i];
        }

        jdoubleArray array = env->NewDoubleArray(size);
        env->SetDoubleArrayRegion(array, 0, size, body);

        return array;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT jlongArray JNICALL Java_xmipp_MetaData_findObjects(JNIEnv *env,
        jobject jobj)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);

    try
    {
        std::vector < size_t > ids;
        md->findObjects(ids);

        // Copies vector into array.
        size_t size = ids.size();
        jlong *body = new jlong[size];
        for (int i = 0; i < size; i++)
        {
            body[i] = ids[i];
        }

        // Sets array value
        jlongArray array = env->NewLongArray(size);
        env->SetLongArrayRegion(array, 0, size, body);

        return array;
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return NULL;
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_importObjects
(JNIEnv *env, jobject jobj, jobject from, jlongArray jids)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    MetaData * mdfrom = GET_INTERNAL_METADATA(from);

    try
    {
        jlong *ids = env->GetLongArrayElements(jids, 0);
        int size = env->GetArrayLength(jids);

        std::vector<size_t> out_ids(size);
        for(int i=0; i<size; i++)
        {
            out_ids[i] = ids[i];
        }

        env->ReleaseLongArrayElements(jids, ids, 0);

        md->importObjects(*mdfrom, out_ids);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jlong JNICALL Java_xmipp_MetaData_firstObject(JNIEnv *env,
        jobject jobj)
{
    std::string msg = "";
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    jlong id = 0;

    try
    {
        id = md->firstObject();
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return id;
}

JNIEXPORT jlong JNICALL Java_xmipp_MetaData_addObject(JNIEnv *env, jobject jobj)
{
    jlong id = 0;
    std::string msg = "";
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    try
    {
        id = md->addObject();
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return id;
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_addLabel(JNIEnv *env, jobject jobj, jint label)
{
    std::string msg = "";
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    try
    {
        md->addLabel((MDLabel)label);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_getPCAbasis
(JNIEnv *env, jobject jmetadata, jobject jbasis)
{
    std::string msg = "";

    try
    {
        MetaData * MDin = GET_INTERNAL_METADATA(jmetadata);
        ImageGeneric *basis= GET_INTERNAL_IMAGE_GENERIC(jbasis);

        MultidimArray<double> *mdArray;
        MULTIDIM_ARRAY_GENERIC(*basis).getMultidimArrayPointer(mdArray);

        ProgAnalyzeCluster program;
        program.NPCA=4;
        program.Niter=10;
        program.dontMask=false;
        program.SFin=*MDin;
        program.produceSideInfo();
        program.pcaAnalyzer.evaluateZScore(program.NPCA, program.Niter);
        program.produceBasis(*mdArray);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_computeFourierStatistics
(JNIEnv *env, jobject jobj, jstring filename)
{
    std::string msg = "";
    MetaData * MDout = GET_INTERNAL_METADATA(jobj);

    try
    {
        const char * fn = env->GetStringUTFChars(filename, false);
        MetaData MDin(fn);
        getFourierStatistics(MDin, 1, *MDout, true, 2);

        env->ReleaseStringUTFChars(filename, fn);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_enableDebug
(JNIEnv *, jobject)
{
    extern int debug;
    debug=1;
}

JNIEXPORT void JNICALL Java_xmipp_MetaData_readPlain
(JNIEnv * env, jobject jobj, jstring jfile, jstring jcolumns)
{
    std::string msg = "";

    try
    {
        const char * nfile = env->GetStringUTFChars(jfile, false);
        const char * ncolumns = env->GetStringUTFChars(jcolumns, false);
        MetaData * metadata = GET_INTERNAL_METADATA(jobj);
        metadata->readPlain(nfile, ncolumns);

        env->ReleaseStringUTFChars(jfile, nfile);
        env->ReleaseStringUTFChars(jcolumns, ncolumns);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}
