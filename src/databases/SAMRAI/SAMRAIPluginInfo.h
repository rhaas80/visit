// ****************************************************************************
//                               SAMRAIPluginInfo.h
// ****************************************************************************

#ifndef SAMRAI_PLUGIN_INFO_H
#define SAMRAI_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;

// ****************************************************************************
//  Class: SAMRAIDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the SAMRAI plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: miller -- generated by xml2info
//  Creation:   Tue Jul 29 08:15:04 PDT 2003
//
//  Modifications:
//
// ****************************************************************************

class SAMRAIGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class SAMRAICommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual SAMRAIGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class SAMRAIMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual SAMRAICommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class SAMRAIEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual SAMRAICommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

#endif
