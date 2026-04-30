-- blake3(input <string>)

local band,bxor,bor,rsh,rrot=bit32.band,bit32.bxor,bit32.bor,bit32.rshift,bit32.rrotate

local IV={
    0x6A09E667,0xBB67AE85,0x3C6EF372,0xA54FF53A,
    0x510E527F,0x9B05688C,0x1F83D9AB,0x5BE0CD19,
}

-- same permutation applied each of the 7 rounds
-- 0-based indices for easier indexing into Lua tables
local PERM={3,7,4,11,8,1,5,14,2,12,13,6,10,15,16,9}

local F_CHUNK_START,F_CHUNK_END,F_PARENT,F_ROOT=0x01,0x02,0x04,0x08
local BLEN,CLEN=64,1024

local function u32(n) return band(n,0xFFFFFFFF) end

local function G(v,a,b,c,d,x,y)
    v[a]=u32(v[a]+v[b]+x); v[d]=rrot(bxor(v[d],v[a]),16)
    v[c]=u32(v[c]+v[d]);   v[b]=rrot(bxor(v[b],v[c]),12)
    v[a]=u32(v[a]+v[b]+y); v[d]=rrot(bxor(v[d],v[a]),8)
    v[c]=u32(v[c]+v[d]);   v[b]=rrot(bxor(v[b],v[c]),7)
end

local function compress(cv,m,counter,blen,flags)
    local v={
        cv[1],cv[2],cv[3],cv[4],
        cv[5],cv[6],cv[7],cv[8],
        IV[1],IV[2],IV[3],IV[4],
        u32(counter),0,blen,flags,
    }
    local w={}; for i=1,16 do w[i]=m[i] or 0 end

    for _=1,7 do
        -- columns
        G(v,1,5, 9,13,w[1], w[2]);  G(v,2,6,10,14,w[3], w[4])
        G(v,3,7,11,15,w[5], w[6]);  G(v,4,8,12,16,w[7], w[8])
        -- diagonals
        G(v,1,6,11,16,w[9], w[10]); G(v,2,7,12,13,w[11],w[12])
        G(v,3,8, 9,14,w[13],w[14]); G(v,4,5,10,15,w[15],w[16])

        local p={}; for i=1,16 do p[i]=w[PERM[i]] end; w=p
    end

    local out={}; for i=1,8 do out[i]=bxor(v[i],v[i+8]) end
    return out
end

-- decode a <64-byte string into 16 little-endian u32 words
local function le_words(s)
    local w={}
    for i=0,15 do
        local o=i*4
        w[i+1]=(string.byte(s,o+1) or 0)
              +((string.byte(s,o+2) or 0)*0x100)
              +((string.byte(s,o+3) or 0)*0x10000)
              +((string.byte(s,o+4) or 0)*0x1000000)
    end
    return w
end

-- compress one chunk (<1024 bytes) down to an 8-word chaining value
local function chunkCV(chunk,iv,counter,rf)
    local n=math.max(1,math.ceil(#chunk/BLEN))
    local cv={table.unpack(iv)}
    for i=1,n do
        local blk=string.sub(chunk,(i-1)*BLEN+1,i*BLEN)
        local f=0
        if i==1 then f=bor(f,F_CHUNK_START) end
        if i==n then f=bor(f,bor(F_CHUNK_END,rf)) end
        cv=compress(cv,le_words(blk),counter,#blk,f)
    end
    return cv
end

-- merge two chaining values as a parent node
local function parentCV(l,r,iv,rf)
    local m={}
    for i=1,8 do m[i]=l[i]; m[i+8]=r[i] end
    return compress(iv,m,0,BLEN,bor(F_PARENT,rf))
end

local function encode_hex(cv)
    local out={}
    for _,word in ipairs(cv) do
        for j=0,3 do out[#out+1]=string.format("%02x",band(rsh(word,j*8),0xFF)) end
    end
    return table.concat(out)
end

local function blake3(input)
    -- split into chunks of CLEN bytes
    local chunks={}
    if #input==0 then
        chunks[1]=""
    else
        for i=1,math.ceil(#input/CLEN) do
            chunks[i]=string.sub(input,(i-1)*CLEN+1,i*CLEN)
        end
    end

    -- compress each chunk; single chunk inputs get ROOT right here
    local single=#chunks==1
    local cvs={}
    for i,chunk in ipairs(chunks) do
        cvs[i]=chunkCV(chunk,IV,i-1,single and F_ROOT or 0)
    end

    -- fold up the merkle tree; ROOT goes on the final merge
    while #cvs>1 do
        local next={}
        for i=1,#cvs,2 do
            if cvs[i+1] then
                next[#next+1]=parentCV(cvs[i],cvs[i+1],IV,#cvs==2 and F_ROOT or 0)
            else
                next[#next+1]=cvs[i] -- odd node promoted unchanged
            end
        end
        cvs=next
    end

    return encode_hex(cvs[1])
end

return blake3
