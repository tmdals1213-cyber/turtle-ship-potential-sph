"""Unchanged generic XML construction helpers extracted from the KCS preparer."""
import xml.etree.ElementTree as ET

def el(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {key: str(value) for key, value in attrs.items()})

def box(parent, point, size):
    b = el(parent, 'drawbox'); el(b, 'boxfill').text = 'solid'
    el(b, 'point', **dict(zip('xyz', point))); el(b, 'size', **dict(zip('xyz', size)))
