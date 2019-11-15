#include <QMouseEvent>
#include <QResizeEvent>
#include <QToolTip>
#include <QPainter>
#include <QMenu>
#include <iostream>
#include "notelabel.hh"
#include "notegraphwidget.hh"

namespace {
	static const int text_margin = 3; // Margin of the label texts
}

const int NoteLabel::render_delay = 150; // How many ms to wait before updating pixmap after some action
const int NoteLabel::resize_margin = 5; // How many pixels is the resize area
const double NoteLabel::default_length = 0.5; // The preferred size of notes
const double NoteLabel::min_length = 0.05; // How many seconds minimum

NoteLabel::NoteLabel(const Note &note, QWidget *parent, bool floating)
	: QLabel(parent), m_note(note), m_selected(false), m_floating(floating), m_resizing(0), m_hotspot()
{
	updateLabel();
	setMouseTracking(true);
	hide();
	// We don't want to show the widget and create the pixmap as that is slow.
	// Since the undo-framework relies on rapidly creating and deleting NoteLabels,
	// this is a necessity to get adequete performance. NoteGraphWidget creates the
	// pixmaps later on once the final NoteLabels have been found.
}

void NoteLabel::updatePixmap()
{
	if (isHidden()) return;
	QFont font;
	font.setStyleStrategy(QFont::ForceOutline);
	QFontMetrics metric(font);
	NoteGraphWidget *ngw = qobject_cast<NoteGraphWidget*>(parent());
	QSize size(100, metric.size(Qt::TextSingleLine, lyric()).height() + 2 * text_margin);
	if (ngw) size.setWidth(ngw->s2px(m_note.length()));

	if (size.isEmpty()) return;

	QImage image(size.width(), size.height(), QImage::Format_ARGB32_Premultiplied);
	image.fill(qRgba(0, 0, 0, 0));

	QLinearGradient gradient(0, 0, 0, image.height()-1);
	float ff = m_floating ? 1.0f : 0.6f;
	int alpha = m_floating ? 160 : ( isSelected() ? 80 : 220 );
	gradient.setColorAt(0.0, m_floating ? QColor(255, 255, 255, alpha) : QColor(50, 50, 50, alpha));
	if (m_note.type == Note::NORMAL) {
		gradient.setColorAt(0.2, QColor(100 * ff, 100 * ff, 255 * ff, alpha));
		gradient.setColorAt(0.8, QColor(100 * ff, 100 * ff, 255 * ff, alpha));
		gradient.setColorAt(1.0, QColor(100 * ff, 100 * ff, 200 * ff, alpha));
	} else if (m_note.type == Note::GOLDEN) {
		gradient.setColorAt(0.2, QColor(255 * ff, 255 * ff, 100 * ff, alpha));
		gradient.setColorAt(0.8, QColor(255 * ff, 255 * ff, 100 * ff, alpha));
		gradient.setColorAt(1.0, QColor(160 * ff, 160 * ff, 100 * ff, alpha));
	} else if (m_note.type == Note::FREESTYLE) {
		gradient.setColorAt(0.2, QColor(100 * ff, 180 * ff, 100 * ff, alpha));
		gradient.setColorAt(0.8, QColor(100 * ff, 180 * ff, 100 * ff, alpha));
		gradient.setColorAt(1.0, QColor(100 * ff, 120 * ff, 100 * ff, alpha));
	}

	{
		QPainter painter(&image);
		painter.setRenderHint(QPainter::Antialiasing);
		painter.setPen(isSelected() ? Qt::red : Qt::black); // Hilight selected note
		painter.setBrush(gradient);
		painter.drawRoundedRect(QRectF(0.5, 0.5, image.width()-1, image.height()-1), 8, 8);

		painter.setFont(font);
		painter.setPen(isSelected() ? Qt::red : Qt::white);
		painter.drawText(QRect(QPoint(text_margin, text_margin), QSize(size.width()-text_margin*2, size.height()-text_margin*2)), Qt::AlignCenter, lyric());

		// Render sentence end indicator
		if (m_note.lineBreak) {
			painter.setPen(QPen(QBrush(QColor(255, 0, 0)), 4));
			painter.drawLine(2, 0, 2, image.height()-1);
		}
	}

	setPixmap(QPixmap::fromImage(image));
	updateTips();
	show();
}

void NoteLabel::setSelected(bool state) {
	if (m_selected != state) {
		m_selected = state;
		QTimer::singleShot(50, this, SLOT(updatePixmap())); // Quick update needed here for box selection
		if (!m_selected) {
			startResizing(0); // Reset
			startDragging(QPoint()); // Reset
		}
	}
}

void NoteLabel::resizeEvent(QResizeEvent *) { updatePixmap(); }

void NoteLabel::moveEvent(QMoveEvent *) { updateTips(); }

void NoteLabel::mouseMoveEvent(QMouseEvent *event)
{
	NoteGraphWidget* ngw = qobject_cast<NoteGraphWidget*>(parent());
	if (m_resizing != 0 && ngw) {
		// Resizing
		double diffsecs = ngw->px2s(event->pos().x());
		if (m_resizing < 0) m_note.begin += diffsecs;  // Left side
		else m_note.end += diffsecs - m_note.length(); // Right side
		// Enforce minimum size
		if (m_note.length() < min_length) {
			if (m_resizing < 0) m_note.begin = m_note.end - min_length; // Left side
			else m_note.end = m_note.begin + min_length; // Right side
		}
		updateLabel();
		ngw->updateNotes(m_resizing > 0);

	} else if (!m_hotspot.isNull() && ngw) {
		// Moving
		QPoint newpos = pos() + event->pos() - m_hotspot;
		double ds = ngw->px2s(event->pos().x() - m_hotspot.x());
		int dn = ngw->px2n(event->pos().y()) - ngw->px2n(m_hotspot.y());
		NoteLabels& labels = ngw->selectedNotes();
		for (int i = 0; i < labels.size(); ++i) {
			NoteLabel *nl = labels[i];
			nl->note().begin += ds;
			nl->note().end += ds;
			nl->note().note += dn;
			nl->updateLabel();
		}
		ngw->updateNotes((event->pos() - m_hotspot).x() < 0);
		// Check if we need a new hotspot, because the note was constrained
		if (pos().x() != newpos.x()) m_hotspot = event->pos();

	} else {
		// Hover cursors
		if (event->pos().x() < NoteLabel::resize_margin || event->pos().x() > width() - NoteLabel::resize_margin) {
			setCursor(QCursor(Qt::SizeHorCursor));
		} else {
			setCursor(QCursor(Qt::OpenHandCursor));
		}
	}
	QToolTip::showText(event->globalPos(), toolTip(), this);
	event->ignore(); // Propagate event to parent
}

void NoteLabel::startResizing(int dir)
{
	m_resizing = dir;
	m_hotspot = QPoint(); // Reset
	if (dir != 0) setCursor(QCursor(Qt::SizeHorCursor));
	else setCursor(QCursor());
}

void NoteLabel::startDragging(const QPoint& point)
{
	m_hotspot = point;
	m_resizing = 0;

	if (!point.isNull()) setCursor(QCursor(Qt::ClosedHandCursor));
	else setCursor(QCursor());
}

void NoteLabel::updateLabel()
{
	NoteGraphWidget* ngw = qobject_cast<NoteGraphWidget*>(parent());
	if (ngw) {
		// Update label geometry
		resize(ngw->s2px(m_note.length()), height());
		move(ngw->s2px(m_note.begin), ngw->n2px(m_note.note) - height() / 2);
	}
}

void NoteLabel::updateTips()
{
	setToolTip(description(true));
	setStatusTip(description(false));
	setWhatsThis(description(true));
}

QString NoteLabel::description(bool multiline) const
{
	MusicalScale ms;
	return QString("Syllable: \"%2\"%1Type: %3%1Note: %4 (%5)%1%6 s - %7 s (= %8 s)")
		.arg(multiline ? "\n" : ", ")
		.arg(lyric())
		.arg(m_note.typeString())
		.arg(ms.getNoteStr(ms.getNoteFreq(m_note.note)))
		.arg(m_note.note)
		.arg(QString::number(m_note.begin, 'f', 4))
		.arg(QString::number(m_note.end, 'f', 4))
		.arg(QString::number(m_note.length(), 'f', 4)
		);
}

NoteLabel::operator Operation() const
{
	Operation op("NEW", -1); // -1 for id means auto-calculate based on position
	op << m_note.syllable << m_note.begin << m_note.end << m_note.note << m_floating << m_note.lineBreak << m_note.getTypeInt();
	return op;
}
