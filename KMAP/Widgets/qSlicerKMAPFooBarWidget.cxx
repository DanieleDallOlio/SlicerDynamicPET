/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerKMAPFooBarWidget.h"
#include "ui_qSlicerKMAPFooBarWidget.h"

//-----------------------------------------------------------------------------
class qSlicerKMAPFooBarWidgetPrivate
  : public Ui_qSlicerKMAPFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerKMAPFooBarWidget);
protected:
  qSlicerKMAPFooBarWidget* const q_ptr;

public:
  qSlicerKMAPFooBarWidgetPrivate(
    qSlicerKMAPFooBarWidget& object);
  virtual void setupUi(qSlicerKMAPFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerKMAPFooBarWidgetPrivate
::qSlicerKMAPFooBarWidgetPrivate(
  qSlicerKMAPFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerKMAPFooBarWidgetPrivate
::setupUi(qSlicerKMAPFooBarWidget* widget)
{
  this->Ui_qSlicerKMAPFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerKMAPFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerKMAPFooBarWidget
::qSlicerKMAPFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerKMAPFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerKMAPFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerKMAPFooBarWidget
::~qSlicerKMAPFooBarWidget()
{
}
