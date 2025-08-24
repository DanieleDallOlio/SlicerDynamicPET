#include "KeyPressWatcher.h"
#include <iostream>

KeyPressWatcher::KeyPressWatcher(QObject* parent)
  : QObject(parent), active_(true)
{
}

void KeyPressWatcher::setActive(bool active)
{
  active_ = active;
}

bool KeyPressWatcher::isActive() const
{
  return active_;
}

bool KeyPressWatcher::eventFilter(QObject* obj, QEvent* event)
{
  if (!active_)
    return QObject::eventFilter(obj, event);


  if (event->type() == QEvent::KeyPress)
  {
    QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
    int key = keyEvent->key();
    if (key == Qt::Key_Delete || key == Qt::Key_Backspace)
    {
      emit deletePressed();
      return true; // consume event
    }
  }
  return QObject::eventFilter(obj, event);
}
