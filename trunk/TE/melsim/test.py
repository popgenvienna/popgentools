#!/usr/bin/env python
import re

s="Dsil\Loa"
famid=re.sub(r"D[^\\]+\\","",s)
print famid
