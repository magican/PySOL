import bufr
f ='/tmp/ascat_20150531_013300_metopa_44686_eps_o_coa_ovw.l2_bufr'
hh = bufr.BUFRFile(f)
var = []
unit = []
data = []
for record in hh:
    for entry in record:
        entryname = (entry.name.lower().replace(' ','_').
                                 strip('_').strip('*'))
        if entryname not in var:
            var.append(entryname)
            unit.append(entry.unit)
            data.append(entry.data[0])
            print entryname,entry.unit,entry.data[0]