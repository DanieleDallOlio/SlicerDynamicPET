/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifdef _WIN32
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif


// KMAP Logic includes
#include <vtkSlicerKMAPLogic.h>

// KMAP includes
#include "qSlicerKMAPModule.h"
#include "qSlicerKMAPModuleWidget.h"



//-----------------------------------------------------------------------------
class qSlicerKMAPModulePrivate
{
public:
  qSlicerKMAPModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerKMAPModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerKMAPModulePrivate::qSlicerKMAPModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerKMAPModule methods

//-----------------------------------------------------------------------------
qSlicerKMAPModule::qSlicerKMAPModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerKMAPModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerKMAPModule::~qSlicerKMAPModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerKMAPModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerKMAPModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerKMAPModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerKMAPModule::icon() const
{
  return QIcon(":/Icons/KMAP.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerKMAPModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerKMAPModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerKMAPModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerKMAPModule
::createWidgetRepresentation()
{
  return new qSlicerKMAPModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerKMAPModule::createLogic()
{
  return vtkSlicerKMAPLogic::New();
}
