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


// DynamicPET Logic includes
#include <vtkSlicerDynamicPETLogic.h>

// DynamicPET includes
#include "qSlicerDynamicPETModule.h"
#include "qSlicerDynamicPETModuleWidget.h"



//-----------------------------------------------------------------------------
class qSlicerDynamicPETModulePrivate
{
public:
  qSlicerDynamicPETModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerDynamicPETModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerDynamicPETModulePrivate::qSlicerDynamicPETModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerDynamicPETModule methods

//-----------------------------------------------------------------------------
qSlicerDynamicPETModule::qSlicerDynamicPETModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerDynamicPETModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerDynamicPETModule::~qSlicerDynamicPETModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerDynamicPETModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerDynamicPETModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDynamicPETModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerDynamicPETModule::icon() const
{
  return QIcon(":/Icons/DynamicPET.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerDynamicPETModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDynamicPETModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerDynamicPETModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerDynamicPETModule
::createWidgetRepresentation()
{
  return new qSlicerDynamicPETModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerDynamicPETModule::createLogic()
{
  return vtkSlicerDynamicPETLogic::New();
}
