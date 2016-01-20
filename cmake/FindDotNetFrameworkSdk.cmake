# Set paths and vars for .NET compilers 
# This is hand-rolled because I had problems with the one from SimpleITK

#
# The following variables are set:
#   CSHARP_DOTNET_FOUND
#   CSHARP_DOTNET_COMPILER_${version} eg. "CSHARP_DOTNET_COMPILER_v4.0.30319"
#   CSHARP_DOTNET_VERSION eg. "v4.0.30319"
#   CSHARP_DOTNET_VERSIONS eg. "v2.0.50727, v3.5, v4.0.30319"
#   DotNetFrameworkSdk_USE_FILE
#
#   CSHARP_PROJECT_BUILDER (xbuild/msbuild)

set(framework_dir "C:/Windows/Microsoft.NET/Framework")

set(CSHARP_DOTNET_VERSION "v4.0.30319")
set(CSHARP_DOTNET_VERSIONS "")
set(CSHARP_DOTNET_COMPILER_${CSHARP_DOTNET_VERSION} "${framework_dir}/${CSHARP_DOTNET_VERSION}/csc.exe")
set(CSHARP_PROJECT_BUILDER "${framework_dir}/${CSHARP_DOTNET_VERSION}/MSBuild.exe")

if(EXISTS ${CSHARP_DOTNET_COMPILER_${CSHARP_DOTNET_VERSION}})
	set(CSHARP_DOTNET_FOUND 1)
else()
	set(CSHARP_DOTNET_FOUND 0)
endif()

# Set USE_FILE
get_filename_component( current_list_path ${CMAKE_CURRENT_LIST_FILE} PATH )
set( DotNetFrameworkSdk_USE_FILE ${current_list_path}/UseDotNetFrameworkSdk.cmake )