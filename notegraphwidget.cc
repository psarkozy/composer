#include <QtGui>
#include <iostream>
#include <cmath>
#include "notelabel.hh"
#include "notegraphwidget.hh"

NoteGraphWidget::NoteGraphWidget(QWidget *parent)
	: QWidget(parent), m_resizing()
{
	setAcceptDrops(true);
}

void NoteGraphWidget::clear()
{
	// Clear NoteLabels
	const QObjectList &childlist = children();
	for (QObjectList::const_iterator it = childlist.begin(); it != childlist.end(); ++it) {
		NoteLabel *child = dynamic_cast<NoteLabel*>(*it);
		if (child) child->close();
	}
	m_notes.clear();
}

void NoteGraphWidget::setLyrics(QString lyrics)
{
	QTextStream ts(&lyrics, QIODevice::ReadOnly);
	const int gap = 10;
	int x = gap, y = 50;

	clear();
	while (!ts.atEnd()) {
		QString word;
		ts >> word;
		if (!word.isEmpty()) {
			NoteLabel *wordLabel = new NoteLabel(word, this, QPoint(x, y));
			x += wordLabel->width() + gap;
		}
	}

	rebuildNoteList();
	// Set first and last to non-floating
	m_notes.front()->disableFloating();
	m_notes.back()->disableFloating();

	m_requiredWidth = x + gap;
	updateWidth();
}

void NoteGraphWidget::updateWidth()
{
	setFixedWidth(m_requiredWidth);
}

void NoteGraphWidget::updateNotes()
{
	FloatingGap gap(0);
	rebuildNoteList();
	// Determine gaps between non-floating notes
	for (NoteLabels::iterator it = m_notes.begin(); it != m_notes.end(); ++it) {
		NoteLabel *child = *it; //dynamic_cast<NoteLabel*>(*it);
		if (!child) continue;
		std::cout << "Child, floating: " << child->isFloating() << std::endl;
		if (child->isFloating()) {
			// Add floating note to gap
			gap.addNote(child);
		} else {
			// Fixed note encountered, handle the gap (divide notes evenly into it)
			std::cout << "Gap end, width: " << gap.width() << ", notes: " << gap.notes.size() << std::endl;
			if (!gap.isEmpty()) {
				gap.end = child->x();
				int step = gap.width() / (gap.notes.size() + 1);
				int x = gap.begin + step * 0.5;
				for (NoteLabels::iterator it2 = gap.notes.begin(); it2 != gap.notes.end(); ++it2, x += step) {
					(*it2)->move(x, (*it2)->y());
				}
			}
			// Start a new gap
			gap = FloatingGap(child->x() + child->width());
			std::cout << "New gap from: " << gap.begin << std::endl;
		}
	}
}

void NoteGraphWidget::rebuildNoteList()
{
	m_notes.clear();
	const QObjectList &childlist = children();
	for (QObjectList::const_iterator it = childlist.begin(); it != childlist.end(); ++it) {
		NoteLabel *child = dynamic_cast<NoteLabel*>(*it);
		if (child) m_notes.push_back(child);
	}
	m_notes.sort(cmpNoteLabelPtr);
}

void NoteGraphWidget::dragEnterEvent(QDragEnterEvent *event)
{
	if (event->mimeData()->hasFormat("application/x-notelabel")) {
		if (children().contains(event->source())) {
			event->setDropAction(Qt::MoveAction);
			event->accept();
		} else {
			event->acceptProposedAction();
		}
	} else {
		event->ignore();
	}
}

void NoteGraphWidget::dragMoveEvent(QDragMoveEvent *event)
{
	if (event->mimeData()->hasFormat("application/x-notelabel")) {
		if (children().contains(event->source())) {
			event->setDropAction(Qt::MoveAction);
			event->accept();
		} else {
			event->acceptProposedAction();
		}
		// Update note positions on-the-fly
		updateNotes();
	} else {
		event->ignore();
	}
}

void NoteGraphWidget::dropEvent(QDropEvent *event)
{
	if (event->mimeData()->hasFormat("application/x-notelabel")) {
		const QMimeData *mime = event->mimeData();
		QByteArray itemData = mime->data("application/x-notelabel");
		QDataStream dataStream(&itemData, QIODevice::ReadOnly);

		QString text;
		QPoint offset;
		dataStream >> text >> offset;
		new NoteLabel(text, this, event->pos() - offset, false);
		updateNotes();

		if (event->source() == this) {
			event->setDropAction(Qt::MoveAction);
			event->accept();
		} else {
			event->acceptProposedAction();
		}
	} else {
		event->ignore();
	}
}

void NoteGraphWidget::mousePressEvent(QMouseEvent *event)
{
	NoteLabel *child = static_cast<NoteLabel*>(childAt(event->pos()));
	if (!child)
		return;

	QPoint hotSpot = event->pos() - child->pos();

	// Left Click
	if (event->button() == Qt::LeftButton) {

		// Determine if it is drag or resize
		const int resizeHandleSize = 5;
		if (hotSpot.x() < resizeHandleSize || hotSpot.x() > child->width() - resizeHandleSize) {
			// Start a resize
			m_resizing = child;
			child->setResizing( (hotSpot.x() < resizeHandleSize) ? -1 : 1 );
			std::cout << "resize" << std::endl;

		} else {
			// Initialize Drag

			QByteArray itemData;
			QDataStream dataStream(&itemData, QIODevice::WriteOnly);
			dataStream << child->getText() << QPoint(hotSpot);

			QMimeData *mimeData = new QMimeData;
			mimeData->setData("application/x-notelabel", itemData);
			mimeData->setText(child->getText());

			QDrag *drag = new QDrag(this);
			drag->setMimeData(mimeData);
			drag->setPixmap(*child->pixmap());
			drag->setHotSpot(hotSpot);

			child->disableFloating();
			child->createPixmap();
			child->hide();

			if (drag->exec(Qt::MoveAction | Qt::CopyAction, Qt::CopyAction) == Qt::MoveAction)
				child->close();
			else
				child->show();
		}

	// Right Click
	} else if (event->button() == Qt::RightButton) {

		// Cut the text in a position proportional to the click point
		float relRatio = float(hotSpot.x()) / child->width();
		int cutpos = int(std::ceil(child->getText().length() * relRatio));
		QString firstst = child->getText().left(cutpos);
		QString secondst = child->getText().right(child->getText().length() - cutpos);

		// Create new labels
		NoteLabel *newLabel1 = new NoteLabel(firstst, this,child->pos());
		new NoteLabel(secondst, this, newLabel1->pos() + QPoint(newLabel1->width(), 0));

		// Delete the old one
		child->close();
	}
}

void NoteGraphWidget::mouseReleaseEvent(QMouseEvent *event)
{
	(void)*event;
	m_resizing->setResizing(0);
	m_resizing = NULL;
}

void NoteGraphWidget::wheelEvent(QWheelEvent *event)
{
	NoteLabel *child = static_cast<NoteLabel*>(childAt(event->pos()));
	if (!child)
		return;

	// Figure out new size and apply it
	int neww = child->size().width() + event->delta() * 0.1;
	child->resize(neww, child->size().height());

	event->accept();
}

void NoteGraphWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
	NoteLabel *child = static_cast<NoteLabel*>(childAt(event->pos()));
	if (!child)
		return;

	// Spawn an input dialog
	bool ok;
	QString text = QInputDialog::getText(this, tr("Edit lyric"),
										  tr("Lyric:"), QLineEdit::Normal,
										  child->getText(), &ok);
	if (ok && !text.isEmpty()) {
		child->setText(text);
		child->createPixmap();
	}

	event->accept();
}


void FloatingGap::addNote(NoteLabel* n)
{
	notes.push_back(n);
	end = n->x();
}
