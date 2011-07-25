# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'GLU copyright and license terms'
__index__     = 'General modules'
__order__     = 5
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See this file or run: glu license'
__revision__  = '$Id$'

license = '''\
GLU SOFTWARE LICENSE

Last amended: 2008-02-05

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health &
Human Services. Funded by NCI under Contract N01-CO-12400.

BioInformed LLC and the U.S. Department of Health & Human Services
(“Copyright Holders”) hereby grant to the public a perpetual, irrevocable,
royalty-free non-exclusive, worldwide license to use, copy, distribute,
display and prepare derivatives of the software in source or binary forms,
with or without modification, subject only to the following conditions:

  o Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  o Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  o Neither the name of BioInformed LLC nor the names of its contributors
    nor the U.S. Department of Health & Human Services may be used to
    endorse or promote products derived from this software without specific
    prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
'''

def main():
  import pydoc
  pydoc.pager(license)
