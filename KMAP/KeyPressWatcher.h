#ifndef KeyPressWatcher_h
#define KeyPressWatcher_h

#include <QObject>
#include <QEvent>
#include <QKeyEvent>

class KeyPressWatcher : public QObject
{
  Q_OBJECT

public:
  explicit KeyPressWatcher(QObject* parent = nullptr);
  void setActive(bool active);
  bool isActive() const;

protected:
  bool eventFilter(QObject* obj, QEvent* event) override;

signals:
  void deletePressed();

private:
  bool active_;
};

#endif // KeyPressWatcher_h
